# maf2con 구현 계획서

## 1. 작업 범위

- 대상 저장소: `https://github.com/chulbioinfo/ReAlignPro`
- 로컬 워크스페이스: `/Users/openhl/Documents/New project 2/ReAlignPro`
- 확인 커밋: `6e53951f4d6133a860136601c592998fa1d85409`
- 목표 기능: `realignpro maf2con`
  - MAF/MAF.GZ 입력을 읽고, reference 좌표계 기준 BED3 constrained interval을 출력한다.
  - target을 전체 종으로 선택했을 때, target group 내 major allele/major sequence similarity가 기본 99% threshold 이상인 reference base들을 constrained base로 판정하고, 인접 base를 region으로 merge한다.

## 2. maf2bed 현재 구조 요약

현재 `maf2bed`는 거의 단일 파일 구조이다.

- CLI 설정/검증: `src/realignpro/maf2bed.py`
  - `Maf2BedConfig`
  - `build_arg_parser()`
  - `_validate_and_build_config()`
- MAF I/O 및 block parsing:
  - `open_text_maybe_gzip()`
  - `list_maf_ids()`
  - `parse_maf_block_lines()`
- 핵심 판정 알고리즘:
  - `matrix2var()`
- multiprocessing pipeline:
  - `worker_proc()`
  - `writer_thread()`
  - `maf2bed_multiprocessing()`
- top-level command 연결:
  - `src/realignpro/cli.py`에서 `maf2bed` subcommand를 lazy import한다.

현재 `matrix2var()`의 판정 로직은 다음과 같다.

1. block 안에 모든 target species가 있어야 처리한다.
2. reference sequence를 왼쪽에서 오른쪽으로 훑되, reference gap column은 좌표가 없으므로 skip한다.
3. 각 alignment column에서 target allele 목록과 non-target allele 목록을 만든다.
4. `set(target_nts)`의 길이가 1이면 target들이 모두 같은 allele을 가진 것으로 본다.
5. 그 allele이 non-target allele set에 없으면 hit로 판정한다.
6. hit가 인접하면 strand를 고려해 BED interval로 merge한다.

중요한 점: `maf2bed` 파일 안에 literal `list(set(...))` 호출은 없지만, 핵심 판정이 `set(target_nts)` 기반이다. 이 방식은 allele의 빈도 정보를 버리기 때문에 "major allele이 몇 %인가"를 계산할 수 없다.

## 3. list(set()) / set 기반 판정의 한계

`set()`은 아래 정보를 잃는다.

- allele별 count
- major allele frequency
- 동률 여부
- gap/N/ambiguous base가 전체 target 수에서 차지하는 비율
- 입력 순서

`maf2con`에서는 constrained 판정이 "모두 완전 동일한가"가 아니라 "major allele similarity가 threshold 이상인가"이므로, `set()` 대신 count-preserving tally를 써야 한다.

권장 방식은 `collections.Counter` 또는 4-base 고정 카운터이다.

```python
from collections import Counter

DNA_BASES = {"A", "C", "G", "T"}

def call_major_base(target_nts, threshold=0.99):
    denom = len(target_nts)
    if denom == 0:
        return None, 0.0

    counts = Counter(nt for nt in target_nts if nt in DNA_BASES)
    if not counts:
        return None, 0.0

    top_count = max(counts.values())
    top_bases = [base for base, count in counts.items() if count == top_count]
    ratio = top_count / float(denom)

    if len(top_bases) != 1:
        return None, ratio
    if ratio < threshold:
        return None, ratio

    return top_bases[0], ratio
```

위 설계에서는 gap, `N`, 기타 ambiguous base는 major 후보에는 넣지 않지만 denominator에는 남긴다. 따라서 missing/gap이 많은 column이 가짜 constrained region으로 판정되는 것을 막는다.

Threshold 비교는 `ratio >= threshold`로 한다. 따라서 기본값 `0.99`에서는 `990/1000 = 99%`도 constrained base로 포함된다.

## 4. maf2con 알고리즘 제안

### 4.1 입력 옵션

새 subcommand:

```bash
realignpro maf2con --input merged.maf.gz --ref-id hg38 --target-ids all --output constrained.bed
```

권장 옵션:

- `--input`: MAF 또는 MAF.GZ
- `--output`: BED3 출력. 생략 시 `*.con.bed`로 derive
- `--threads`: `maf2bed`와 동일하게 reader/writer 포함 total threads
- `--ref-id`: BED 좌표계 기준 species ID
- `--target-ids`: comma-separated list 또는 `all`
- `--outgroup-ids`: constrained 계산에서 제외할 species
- `--min-major-similarity`: 기본 `0.99`
- `--min-target-count`: 기본 `2`
- `--start-method`, `--work-qsize`, `--out-qsize`: `maf2bed`와 동일

`--target-ids all` 처리 방식:

1. 실행 시작 시 `list_maf_ids(input_maf)`로 전체 species ID를 pre-scan한다.
2. `outgroup_ids`를 제외한다.
3. 남은 ID 전체를 target group으로 확정한다.
4. 모든 block은 기존 `maf2bed`처럼 target group 전체가 block 안에 있을 때만 처리한다.

이 방식은 "전체 종 기준 constrained"를 block마다 다른 종 수로 느슨하게 판정하는 문제를 피한다. 단점은 `all` 모드에서 입력을 한 번 더 읽는 비용이 있다는 점이다. 큰 gz MAF에서는 비용이 있지만, constrained 판정의 정확성과 재현성을 우선하면 이 선택이 안전하다.

### 4.2 block별 constrained base 판정

새 핵심 함수는 `matrix2con()`으로 둔다.

입력:

- `block`
- `ref_id`
- `target_ids`
- `outgroup_ids`
- `min_major_similarity`
- `min_target_count`

처리:

1. `target_ids` 전체가 block에 없으면 skip한다.
2. `len(target_ids) < min_target_count`이면 skip한다.
3. reference sequence를 훑으며 reference gap column은 skip한다.
4. 각 reference base column에서 target base를 `target_ids` 순서대로 수집한다.
5. missing 또는 out-of-range sequence는 `"-"`로 처리한다.
6. `A/C/G/T`만 major 후보로 count한다.
7. denominator는 전체 target species 수로 둔다.
8. top allele이 하나이고, `top_count / len(target_ids) >= min_major_similarity`이면 constrained base로 판정한다.
9. constrained base들은 기존 `matrix2var()`의 strand-aware merge 로직을 재사용한다.
10. output은 BED3 line으로 쓴다.

`maf2bed`와 달리 non-target allele absence 조건은 쓰지 않는다. target이 전체 종이면 non-target group이 없고, `maf2con`의 목적도 "target 내부 보존성"이기 때문이다.

### 4.3 Pseudocode

```python
def matrix2con(block, ref_id, target_ids, outgroup_ids, min_major_similarity, min_target_count):
    targets = [sid for sid in target_ids if sid not in set(outgroup_ids)]

    if len(targets) < min_target_count:
        return []
    if not set(targets).issubset(block.keys()):
        return []

    ref = block[ref_id]
    init_reference_coordinate_cursor(ref)

    for aln_idx, ref_nt in enumerate(ref["seq"]):
        if ref_nt == "-":
            continue

        advance_reference_coordinate()

        target_nts = []
        for sid in targets:
            seq = block[sid]["seq"]
            nt = seq[aln_idx].upper() if aln_idx < len(seq) else "-"
            target_nts.append(nt)

        major_base, ratio = call_major_base(target_nts, min_major_similarity)
        is_constrained = major_base is not None

        update_current_interval_or_flush(is_constrained)

    flush_current_interval()
    return bed_lines
```

## 5. 구현 파일 계획

최소 변경 경로:

1. `src/realignpro/maf2con.py` 신규 생성
   - `maf2bed.py`에서 안정적으로 재사용 가능한 helper를 import한다.
   - 예: `open_text_maybe_gzip`, `list_maf_ids`, `parse_maf_block_lines`, queue constants 또는 자체 constants
   - `Maf2ConConfig`, parser, validator, `call_major_base()`, `matrix2con()`, worker/writer/orchestrator/main을 둔다.
2. `src/realignpro/cli.py` 수정
   - help text에 `maf2con` 추가
   - `if cmd == "maf2con": from . import maf2con`
3. `src/realignpro/__init__.py` 문구 업데이트
4. `README.md` quick start 및 기능 목록 업데이트
5. `CHANGELOG.md`에 신규 기능 예정/추가 항목 작성
6. `tests/test_help.py`에 `maf2con --help` smoke test 추가
7. `tests/test_maf2con.py` 신규 추가
   - `call_major_base()` unit tests
   - `matrix2con()` unit tests

나중에 정리할 수 있는 경로:

- `maf2bed.py`와 `maf2con.py`가 중복되는 worker/writer/orchestration 코드를 갖게 되면, 2차 리팩터링으로 `src/realignpro/maf_common.py`를 만들 수 있다.
- 첫 구현에서는 리팩터링 범위를 줄이기 위해 신규 파일 중심으로 가는 편이 안전하다.

## 6. 테스트 계획

### 6.1 단위 테스트

`call_major_base()`:

- all A: constrained
- 99% 이상/미만 threshold 경계
- A/C 동률: not constrained
- gap 또는 N 포함: denominator penalty로 비보존 처리 가능
- empty target list: not constrained

`matrix2con()`:

- plus strand에서 인접 constrained base가 하나의 interval로 merge되는지 확인
- plus strand에서 중간 mismatch가 있으면 interval이 나뉘는지 확인
- minus strand에서 좌표가 forward-strand BED interval로 merge되는지 확인
- target 하나가 block에 없으면 skip되는지 확인
- `outgroup_ids`가 target에서 제외되는지 확인

### 6.2 CLI 테스트

- `python -m realignpro maf2con --help`
- 작은 MAF fixture로 `realignpro maf2con` 실행 후 BED3 출력 비교
- `--target-ids all`이 전체 ID pre-scan 후 outgroup 제외를 적용하는지 확인

### 6.3 회귀 테스트

- 기존 `maf2bed`, `fa2maf`, `tsv2fig` help smoke test가 깨지지 않는지 확인
- `pytest` 전체 실행

## 7. 주요 결정 사항

1. Threshold 비교
   - `ratio >= threshold`
   - 기본값 `0.99`에서는 `990/1000 = 99%`도 constrained로 포함한다.

2. Gap/N 처리
   - 권장: major 후보에서는 제외하지만 denominator에는 포함한다.
   - 이유: gap/N이 많은 column이 constrained로 부풀려지는 것을 방지한다.

3. `--target-ids all`의 의미
   - 권장: 파일 전체에서 발견한 모든 species에서 outgroup을 뺀 목록.
   - block-local all은 block마다 denominator가 달라져 결과 해석이 흔들릴 수 있으므로 피한다.

4. 출력 형식
   - 1차 구현은 `maf2bed`와 동일하게 BED3.
   - 추후 필요하면 BED4/BED5로 major ratio 또는 target count를 score/name에 담는 옵션을 추가할 수 있다.

## 8. 예상 리스크

- `--target-ids all`은 input을 pre-scan하므로 큰 gz MAF에서 시간이 추가된다.
- 현재 MAF parser는 같은 species ID가 한 block에 여러 번 나오면 마지막 record가 이전 것을 덮어쓴다. 이는 기존 `maf2bed`와 같은 동작이므로 이번 기능의 직접 범위 밖이지만, multi-chain MAF를 다룰 때는 주의가 필요하다.
- 99% 이상 조건은 species 수가 작으면 사실상 100% 일치와 같아진다. 예를 들어 target 10개에서 9/10은 90%, 10/10만 100%다.
- target 전체가 모든 block에 있어야 하므로, 일부 species가 빠지는 block은 constrained로 출력되지 않는다. "block에 존재하는 species만으로 판정"이 필요하면 별도 옵션이 필요하다.

## 9. 권장 구현 순서

1. `maf2con.py` 신규 파일 작성
2. `call_major_base()`와 `matrix2con()` 단위 테스트 작성
3. CLI 연결 및 help smoke test 추가
4. 작은 MAF fixture로 end-to-end 테스트 작성
5. README/CHANGELOG 문서 업데이트
6. `pytest`로 검증
