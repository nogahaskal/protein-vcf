import json
import click
import re
from pydantic import BaseModel
from typing import Dict, List, Optional

LETTER_REGEX = r"^(\d+)(.*)"
INSERT_OPERATOR = ">"
DEL_OPERATOR = "del"
DEL_STRING = "<DEL>"


class Diff(BaseModel):
    pos: str
    ref: str
    alt: Optional[str]
    info: Optional[str]


class Haplotype(BaseModel):
    frequency: float
    samples: Dict[str, int]
    aligned_sequences: List[str]
    diffs: List[Diff]
    name: str
    id: str


class Item(BaseModel):
    prot: str
    pos: str
    id: str
    ref: str
    alt: List[str]
    info: str
    samples: Dict[str, List[int]]


HEADER_TO_KEY_MAP = {
    'PROT': 'prot',
    'POS': 'pos',
    'ID': 'id',
    'REF': 'ref',
    'ALT': 'alt',
    'INFO': 'info'
}

ALT_SEPARATOR = ','
SAMPLE_SEPARATOR = '|'

def get_formatted_haplotypes(json_content: dict) -> tuple[List[Haplotype], List[str]]:
    sample_ids = set()
    protein_haplotypes = json_content.get("protein_haplotypes", [])
    formatted_haplotypes = []
    for haplotype in protein_haplotypes:
        frequency = haplotype.get("frequency", 0)
        samples = haplotype.get("samples", {})
        sample_ids.update(samples.keys())

        id = haplotype.get("hex", "")
        diffs = []
        aligned_sequences = haplotype.get("aligned_sequences", [""])
        for diff_data in haplotype.get("diffs", []):
            diff = diff_data.get("diff")
            formatted_diff = format_diff(diff, aligned_sequences[0])
            if not formatted_diff:
                continue
            diffs.append(formatted_diff)

        name = haplotype.get("name", "")
        formatted_haplotype = Haplotype(
            frequency=frequency,
            samples=samples,
            aligned_sequences=aligned_sequences,
            diffs=diffs,
            name=name,
            id=id,
        )
        formatted_haplotypes.append(formatted_haplotype)

    return formatted_haplotypes, list(sample_ids)


def format_diff(diff: str | None, ref_sequence: str) -> Diff | None:
    if not diff:
        return None

    matches = re.match(LETTER_REGEX, diff)
    if matches:
        pos = matches.group(1)
        remaining_part = matches.group(2)
        if INSERT_OPERATOR in remaining_part:
            [ref, alt] = remaining_part.split(INSERT_OPERATOR)
            return Diff(pos=pos, ref=ref, alt=alt, info="")

        if DEL_OPERATOR in remaining_part:
            [_, deletion_part] = remaining_part.split(DEL_OPERATOR)
            ref = ref_sequence[int(pos)]
            alt = DEL_STRING
            info = (
                f"SVTYPE=DEL;END={int(pos) + int(deletion_part.strip('{').strip('}'))}"
            )
            return Diff(pos=pos, ref=ref, alt=alt, info=info)

        return None
    return None


def handle_samples(prev_samples: dict, haplotype: Haplotype, alt_index: int):
    samples: Dict[str, List] = {**prev_samples}

    for sample_id, value in haplotype.samples.items():
        prev_sample = prev_samples.get(sample_id)
        samples_value = [prev_sample] if prev_sample else []
        for _ in range(value):
            samples_value.append(alt_index)
            samples[sample_id] = samples_value

    return samples


def build_items(formatted_haplotypes: List[Haplotype]) -> List[Item]:
    vcf_rows = {}
    for haplotype in formatted_haplotypes:
        for diff in haplotype.diffs:
            row_key = f"{diff.pos}:{diff.ref}"
            prev_row = dict(vcf_rows.get(row_key, {}))

            prot = haplotype.name.split(":")[0]
            pos = diff.pos
            ref = diff.ref
            info = diff.info
            id = haplotype.id
            current_alt = diff.alt

            alt: List[str] = prev_row.get("alt", [])
            if current_alt not in alt:
                alt.append(current_alt)

            prev_samples = prev_row.get("samples", {})
            samples: Dict[str, List] = handle_samples(
                prev_samples, haplotype, alt_index=len(alt)
            )

            vcf_rows[row_key] = Item(
                pos=pos,
                prot=prot,
                id=id,
                ref=ref,
                alt=alt,
                info=info,
                samples=samples,
            )

    return vcf_rows.values()


def append_samples_to_items(
    rows: List[Item], all_sample_ids: List[str]
) -> List[Item]:
    rows_with_samples = []
    for row in rows:
        row.samples = {
            **{sample_id: [0, 0] for sample_id in all_sample_ids},
            **row.samples,
        }
        for sample_id, value in row.samples.items():
            if len(value) < 2:
                row.samples[sample_id].extend([0] * (2 - len(value)))
        rows_with_samples.append(row)
    return rows_with_samples

def generate_vcf_rows(items: List[Item], headers: List[str]):
    rows = []
    for item in items:
        row = []
        dict_item = dict(item)
        for header in headers:
            key = HEADER_TO_KEY_MAP.get(header)
            if key:
                cell = dict_item[key]
                if header == 'ALT':
                    cell = ALT_SEPARATOR.join(cell)
            else:
                sample = item.samples.get(header)
                cell = SAMPLE_SEPARATOR.join(str(num) for num in sample)
            
            row.append(cell)
        rows.append(row)
    return rows
        



def build_vcf_file(vcf_rows: List[List[str]], output: str):
    for row in vcf_rows:
        print(row)


@click.group("pvcf")
def pvcf():
    pass


@pvcf.command("convert")
@click.option("--path", "path", "-p", required=True)
@click.option("--output", "output", "-o", required=True)
def convert_json_to_pvcf(path: str, output: str):
    with open(path, "r") as file:
        try:
            json_content = json.load(file)
            formatted_haplotypes, sample_ids = get_formatted_haplotypes(json_content)
            items = build_items(formatted_haplotypes)
            items_with_samples = append_samples_to_items(items, sample_ids)
            vcf_headers = [*HEADER_TO_KEY_MAP.keys(), *[str(sample_id) for sample_id in sample_ids]]
            vcf_rows = generate_vcf_rows(items_with_samples, vcf_headers)

        except Exception as e:
            print(e)


if __name__ == "__main__":
    pvcf()
