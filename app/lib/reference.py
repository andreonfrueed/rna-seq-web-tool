"""参考文件（基因组 FASTA + GTF）下载与校验。

下载策略：FASTA 用 Ensembl `current_fasta` 稳定文件名（primary_assembly，
比 toplevel 少了补丁/单倍型序列，RNA-seq 比对更干净）；
GTF 文件名含版本号，先抓 `current_gtf` 目录列表解析最新文件名，
失败则退回 web_config 里的固定版本号。

完整性校验（v2.1 新增）：下载后与 Ensembl 同目录 CHECKSUMS 文件交叉校验
（BSD sum 校验和 + 文件大小）。校验不通过直接报错并重下，
防止损坏或被替换的参考文件悄悄参与比对导致整批结果错误。
CHECKSUMS 拿不到时（镜像站没有该文件等）跳过校验，仅靠格式校验兜底。

复用策略：只要物种目录里已有合法的 `.fa` 和 `.gtf`（不管文件名），直接复用，
这样"网页下载过"和"用户手动放文件"都适用。校验是真正的内容校验：
FASTA 必须以 '>' 开头且足够大，GTF 必须有 gene_id 属性，
两者染色体命名风格（chr1 vs 1）必须一致——下了一半的文件不会被当成好的复用。
目录里有多个候选文件时不再随机取用（glob 顺序不确定，可能选中陈旧/错误的
基因组）：FASTA 优先选与官方文件名一致的，GTF 取版本号最高的，
其余情况明确报错要求用户清理。
"""
from __future__ import annotations
import gzip
import re
import shutil
import urllib.request
from pathlib import Path

from .config import load_config

# species -> {fa 子路径（拼在 current_fasta/ 下）, gtf 物种路径段}
REF_SPECS = {
    "hsapiens": {
        "fa_rel": "homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz",
        "gtf_rel_dir": "homo_sapiens",
        "gtf_prefix": "Homo_sapiens.GRCh38.",
    },
    "mmusculus": {
        "fa_rel": "mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz",
        "gtf_rel_dir": "mus_musculus",
        "gtf_prefix": "Mus_musculus.GRCm39.",
    },
}

# 解压后的基因组 FASTA 至少应有这么大（人/小鼠都在 2GB 以上，留足余量）
_MIN_FA_SIZE = 50_000_000


def validate_reference(genome: Path, gtf: Path) -> list[str]:
    """返回错误列表，空列表表示可用。"""
    errors: list[str] = []
    genome, gtf = Path(genome), Path(gtf)
    if not genome.exists() or genome.stat().st_size < _MIN_FA_SIZE:
        errors.append(f"基因组 FASTA 缺失或过小（可能下载不完整）: {genome}")
    if not gtf.exists() or gtf.stat().st_size == 0:
        errors.append(f"注释 GTF 缺失或为空: {gtf}")
    if errors:
        return errors

    # FASTA 内容检查：必须以 '>' 开头
    with open(genome, "rb") as f:
        head = f.read(65536)
    if not head.startswith(b">"):
        errors.append(f"基因组文件不是有效的 FASTA 格式（首字符不是 '>'）: {genome}")
        return errors
    fa_first_seq = head.split(b"\n", 1)[0][1:].split()[0].decode(errors="replace")

    # GTF 内容检查：第一条记录要有 9 列且含 gene_id
    gtf_first_seq = None
    with open(gtf, encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or "gene_id" not in parts[8]:
                errors.append(f"GTF 格式异常（不是 9 列或缺少 gene_id 属性）: {gtf}")
                return errors
            gtf_first_seq = parts[0]
            break
    if gtf_first_seq is None:
        errors.append(f"GTF 文件为空: {gtf}")
        return errors

    # 染色体命名风格必须一致（chr1 vs 1），不一致数基因会全是 0
    if fa_first_seq and gtf_first_seq:
        fa_chr = fa_first_seq.startswith("chr")
        gtf_chr = gtf_first_seq.startswith("chr")
        if fa_chr != gtf_chr:
            errors.append(
                f"基因组与 GTF 的染色体命名风格不一致"
                f"（基因组叫「{fa_first_seq}」，GTF 叫「{gtf_first_seq}」），"
                "混用会导致数不出任何基因，请换成同一来源（建议都用 Ensembl）。")
    return errors


def _expected_fa_name(species: str) -> str:
    return REF_SPECS[species]["fa_rel"].split("/")[-1].replace(".gz", "")


def _genome_file_name(species: str, ref_dir: Path) -> Path:
    return ref_dir / species / _expected_fa_name(species)


def _gtf_file_name(species: str, ref_dir: Path, rel: int) -> Path:
    prefix = REF_SPECS[species]["gtf_prefix"]
    return ref_dir / species / f"{prefix}{rel}.gtf"


def _find_existing(dest_dir: Path, species: str) -> tuple[Path | None, Path | None]:
    """找物种目录里现有的 .fa 与 .gtf（排除 .gz）。

    各只有一个时直接用；多于一个时按确定性规则挑选，绝不随机取用：
    - FASTA：文件名与 REF_SPECS 官方名一致的那个；否则报错要求清理。
    - GTF：符合版本前缀的文件里取 release 号最高的；全不符合则报错。
    """
    fas = sorted(dest_dir.glob("*.fa"))
    gtfs = sorted(dest_dir.glob("*.gtf"))

    fa: Path | None = None
    if len(fas) == 1:
        fa = fas[0]
    elif len(fas) > 1:
        expected = _expected_fa_name(species)
        matches = [p for p in fas if p.name == expected]
        if len(matches) == 1:
            fa = matches[0]
        else:
            names = "、".join(p.name for p in fas)
            raise RuntimeError(
                f"参考目录里找到 {len(fas)} 个基因组文件（{names}），无法确定该用哪个。"
                f"请到 {dest_dir} 里只保留正确的一个 .fa 后重试。")

    gtf: Path | None = None
    if len(gtfs) == 1:
        gtf = gtfs[0]
    elif len(gtfs) > 1:
        prefix = REF_SPECS[species]["gtf_prefix"]
        pat = re.compile(re.escape(prefix) + r"(\d+)\.gtf$")
        best_rel, best = -1, None
        for p in gtfs:
            m = pat.search(p.name)
            if m and int(m.group(1)) > best_rel:
                best_rel, best = int(m.group(1)), p
        if best is None:
            names = "、".join(p.name for p in gtfs)
            raise RuntimeError(
                f"参考目录里找到 {len(gtfs)} 个注释文件（{names}），且都不匹配「{prefix}版本号.gtf」"
                f"命名。请到 {dest_dir} 里只保留正确的一个 .gtf 后重试。")
        gtf = best
    return fa, gtf


def discover_gtf_url(species: str, base: str, listing_text: str) -> str | None:
    prefix = REF_SPECS[species]["gtf_prefix"]
    m = re.search(rf'href="{re.escape(prefix)}(\d+)\.gtf\.gz"', listing_text)
    if not m:
        return None
    rel = m.group(1)
    sub = REF_SPECS[species]["gtf_rel_dir"]
    return f"{base}/current_gtf/{sub}/{prefix}{rel}.gtf.gz"


# ---------------------------------------------------------------- 完整性校验
def _fetch_checksums(url: str) -> dict[str, tuple[int, int]]:
    """解析 Ensembl CHECKSUMS 文件，返回 {文件名: (BSD校验和, 1K块数)}。

    CHECKSUMS 每行格式为 `sum` 命令输出：<校验和> <1K块数> <文件名>。
    拿不到（网络失败/镜像无此文件）返回空字典，调用方跳过校验。
    """
    out: dict[str, tuple[int, int]] = {}
    try:
        text = _http_text(url)
    except Exception:
        return out
    for line in text.splitlines():
        parts = line.split()
        if len(parts) >= 3 and parts[0].isdigit() and parts[1].isdigit():
            out[parts[-1]] = (int(parts[0]), int(parts[1]))
    return out


def _bsd_sum(path: Path) -> tuple[int, int]:
    """BSD sum 算法（与 Ensembl CHECKSUMS 使用的 `sum` 命令一致）。

    返回 (校验和, 文件字节数)。对 ~900MB 的基因组压缩包约需 1-2 分钟，
    仅在首次下载时执行一次。
    """
    r = 0
    total = 0
    with open(path, "rb") as f:
        while True:
            chunk = f.read(1 << 20)
            if not chunk:
                break
            total += len(chunk)
            for b in chunk:
                r = ((r >> 1) + ((r & 1) << 15) + b) & 0xFFFF
    return r, total


def _verify_checksum(gz: Path, checksums: dict[str, tuple[int, int]]) -> None:
    """按 CHECKSUMS 校验已下载文件；不通过则删除文件并抛错（下次会重下）。

    checksums 为空或不含该文件时跳过（无法校验，靠后续格式校验兜底）。
    """
    if not checksums:
        return
    entry = checksums.get(Path(gz).name)
    if entry is None:
        return
    expected_sum, expected_blocks = entry
    actual_sum, actual_bytes = _bsd_sum(Path(gz))
    actual_blocks = (actual_bytes + 1023) // 1024
    if actual_sum != expected_sum or actual_blocks != expected_blocks:
        Path(gz).unlink(missing_ok=True)
        raise RuntimeError(
            f"{Path(gz).name} 完整性校验失败（与 Ensembl CHECKSUMS 不一致："
            f"期望校验和 {expected_sum}，实际 {actual_sum}）。"
            "文件可能损坏或被篡改，已删除，请重新下载。")


def ensure_reference(species: str, ref_dir: Path, progress_cb=None) -> dict[str, str]:
    """确保参考文件就绪，返回 {'genome': str, 'gtf': str}。"""
    if species not in REF_SPECS:
        raise ValueError(f"不支持的物种: {species}")
    ref_dir = Path(ref_dir)
    spec = REF_SPECS[species]
    dest_dir = ref_dir / species
    dest_dir.mkdir(parents=True, exist_ok=True)

    def cb(p, msg):
        if progress_cb:
            progress_cb(min(max(p, 0.0), 1.0), msg)

    existing_fa, existing_gtf = _find_existing(dest_dir, species)
    if existing_fa and existing_gtf:
        if not validate_reference(existing_fa, existing_gtf):
            return {"genome": str(existing_fa), "gtf": str(existing_gtf)}
        cb(0.01, "发现损坏的旧参考文件，已清理")
        existing_fa.unlink(missing_ok=True)
        existing_gtf.unlink(missing_ok=True)
    else:
        # 只找到半边也清掉，避免下载后新旧文件混杂
        for p in (existing_fa, existing_gtf):
            if p:
                p.unlink(missing_ok=True)

    cfg = load_config()
    base = cfg["ensembl_base"].rstrip("/")

    cb(0.02, "开始下载参考文件")

    # FASTA（稳定文件名，用 current_fasta 别名）
    fa_dir_rel = "/".join(spec["fa_rel"].split("/")[:-1])
    fa_url = f"{base}/current_fasta/{spec['fa_rel']}"
    fa_gz = dest_dir / (spec["fa_rel"].split("/")[-1])
    _download(fa_url, fa_gz,
              cb=lambda frac: cb(0.02 + 0.50 * frac, f"下载基因组 FASTA… {int(frac * 100)}%"))
    cb(0.53, "校验基因组完整性（对照 Ensembl CHECKSUMS，约 1-2 分钟）…")
    _verify_checksum(fa_gz, _fetch_checksums(f"{base}/current_fasta/{fa_dir_rel}/CHECKSUMS"))
    genome_path = _genome_file_name(species, ref_dir)
    cb(0.56, "解压基因组 FASTA…")
    _gunzip(fa_gz, genome_path)
    cb(0.60, "基因组 FASTA 就绪")

    # GTF（文件名含版本号 → 先列目录再下，失败退回固定版本）
    gtf_gz = None
    gtf_cks_url = None
    try:
        listing = _http_text(f"{base}/current_gtf/{spec['gtf_rel_dir']}/")
        gtf_url = discover_gtf_url(species, base, listing)
        if gtf_url:
            fname = gtf_url.rsplit("/", 1)[-1]
            gtf_gz = dest_dir / fname
            gtf_cks_url = f"{base}/current_gtf/{spec['gtf_rel_dir']}/CHECKSUMS"
            _download(gtf_url, gtf_gz,
                      cb=lambda frac: cb(0.60 + 0.24 * frac, f"下载注释 GTF… {int(frac * 100)}%"))
    except Exception:
        gtf_gz = None
    if gtf_gz is None or not gtf_gz.exists():
        rel = cfg["ensembl_release"]
        gtf_url = f"{base}/release-{rel}/gtf/{spec['gtf_rel_dir']}/{spec['gtf_prefix']}{rel}.gtf.gz"
        gtf_gz = dest_dir / f"{spec['gtf_prefix']}{rel}.gtf.gz"
        gtf_cks_url = f"{base}/release-{rel}/gtf/{spec['gtf_rel_dir']}/CHECKSUMS"
        _download(gtf_url, gtf_gz,
                  cb=lambda frac: cb(0.60 + 0.24 * frac, f"下载注释 GTF… {int(frac * 100)}%"))

    cb(0.86, "校验注释文件完整性（对照 Ensembl CHECKSUMS）…")
    _verify_checksum(gtf_gz, _fetch_checksums(gtf_cks_url))

    gtf_path = _gtf_file_name(species, ref_dir, _release_of(gtf_gz.name, spec["gtf_prefix"]))
    cb(0.90, "解压注释 GTF…")
    _gunzip(gtf_gz, gtf_path)

    cb(0.96, "校验参考文件…")
    errs = validate_reference(genome_path, gtf_path)
    if errs:
        raise RuntimeError("参考文件校验失败: " + "; ".join(errs))
    cb(1.0, "参考文件准备完成")
    return {"genome": str(genome_path), "gtf": str(gtf_path)}


def _release_of(fname: str, prefix: str) -> int:
    m = re.search(re.escape(prefix) + r"(\d+)\.gtf", fname)
    return int(m.group(1)) if m else load_config()["ensembl_release"]


def _http_text(url: str) -> str:
    with urllib.request.urlopen(url, timeout=60) as r:
        return r.read().decode("utf-8", errors="replace")


def _download(url: str, dest: Path, cb=None) -> None:
    """分块下载并汇报进度；失败时清理半成品，避免下次被当成好文件复用。"""
    dest = Path(dest)
    tmp = dest.with_suffix(dest.suffix + ".part")
    try:
        with urllib.request.urlopen(url, timeout=120) as r, open(tmp, "wb") as f:
            total = int(r.headers.get("Content-Length") or 0)
            got = 0
            while True:
                chunk = r.read(4 * 1024 * 1024)
                if not chunk:
                    break
                f.write(chunk)
                got += len(chunk)
                if cb and total:
                    cb(got / total)
        tmp.replace(dest)
    except Exception:
        tmp.unlink(missing_ok=True)
        dest.unlink(missing_ok=True)
        raise


def _gunzip(src: Path, dest: Path) -> None:
    src, dest = Path(src), Path(dest)
    if not src.name.endswith(".gz"):
        shutil.copyfile(src, dest)
        return
    try:
        with gzip.open(src, "rb") as fin, open(dest, "wb") as fout:
            shutil.copyfileobj(fin, fout)
    except Exception:
        dest.unlink(missing_ok=True)  # 解压到一半的文件不能留
        raise
    src.unlink(missing_ok=True)
