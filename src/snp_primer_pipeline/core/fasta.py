from __future__ import annotations

import logging
from collections.abc import Iterator
from pathlib import Path

logger = logging.getLogger(__name__)


def parse_fasta(file_path: Path) -> Iterator[tuple[str, str]]:
    """Parse a FASTA file and yield (sequence_id, sequence) tuples.

    Args:
        file_path: Path to the FASTA file

    Yields:
        Tuples of (sequence_id, sequence_string)

    Raises:
        OSError: If the file cannot be read
    """
    try:
        with open(file_path) as f:
            current_id: str | None = None
            current_seq: list[str] = []

            for raw_line in f:
                line = raw_line.strip()
                if not line:
                    continue

                if line.startswith(">"):
                    if current_id is not None:
                        yield (current_id, "".join(current_seq))
                    current_id = line[1:]
                    current_seq = []
                else:
                    current_seq.append(line)

            if current_id is not None:
                yield (current_id, "".join(current_seq))

    except OSError:
        logger.exception("Failed to read FASTA file %s", file_path)
        raise


def parse_fasta_to_dict(file_path: Path) -> dict[str, str]:
    """Parse a FASTA file into a dictionary mapping sequence IDs to sequences.

    Args:
        file_path: Path to the FASTA file

    Returns:
        Dictionary of {sequence_id: sequence_string}

    Raises:
        OSError: If the file cannot be read
    """
    return dict(parse_fasta(file_path))
