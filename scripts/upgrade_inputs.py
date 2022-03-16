#!/usr/bin/env python3

import argparse
import dataclasses
import difflib
from typing import List, Optional
import re
import io

import f90nml


@dataclasses.dataclass(unsafe_hash=True)
class NamelistKey:
    namelist: str
    key: str


@dataclasses.dataclass(unsafe_hash=True)
class Replacement:
    old: NamelistKey
    new: NamelistKey


# A list of old namelist options and their new names
REPLACEMENTS = [
    Replacement(
        old=NamelistKey("dist_fn_knobs", "adiabatic_option"),
        new=NamelistKey("physics_flags", "adiabatic_option"),
    )
]


def existing_old_values(filename: str, replacements: dict) -> Optional[dict]:
    """Return any old values present in filename"""
    original = f90nml.read(filename)

    old_values = {}
    for replacement in replacements:
        # Check if old is in original
        namelist = original.get(replacement.old.namelist, False)
        if not namelist:
            continue

        value = namelist.get(replacement.old.key, None)
        if value is not None:
            old_values[replacement] = value

    return old_values


def remove_old_values(filename: str, old_values: dict) -> str:
    """Return the contents of filename with old values removed"""
    old_keys = "|".join((old_value.old.key for old_value in old_values))
    skip_line_re = re.compile(f" *({old_keys}) *=.*")

    new_file = []
    with open(filename, "r") as f:
        for line in f:
            # If this line uses one of the old keys, don't include it in the new file
            if skip_line_re.match(line):
                continue
            new_file.append(line)

    return "".join(new_file)


def apply_fixes(filename: str, replacements: dict):
    """Return a modified version of the input file"""
    old_values = existing_old_values(filename, replacements)
    if old_values == {}:
        return None

    new_file = io.StringIO(remove_old_values(filename, old_values))

    new_values = {
        old_key.new.namelist: {old_key.new.key: value}
        for old_key, value in old_values.items()
    }

    patched_file = io.StringIO()

    f90nml.patch(new_file, new_values, patched_file)
    return patched_file.getvalue()


def yes_or_no(question):
    """Convert user input from yes/no variations to True/False"""
    while True:
        reply = input(question + " [y/N] ").lower().strip()
        if not reply or reply[0] == "n":
            return False
        if reply[0] == "y":
            return True


def create_patch(filename: str, modified: str):
    """Create a unified diff between original and modified"""

    with open(filename, "r") as f:
        original = f.read()

    patch = "\n".join(
        difflib.unified_diff(
            original.splitlines(),
            modified.splitlines(),
            fromfile=filename,
            tofile=filename,
            lineterm="",
        )
    )

    return patch


def possibly_apply_patch(
    filename: str, patch: str, modified: str, quiet=False, force=False
):
    """Possibly apply patch to options_file. If force is True, applies the
    patch without asking, overwriting any existing file. Otherwise,
    ask for confirmation from stdin

    """
    if not quiet:
        print("\n******************************************")
        print("Changes to {}\n".format(filename))
        print(patch)
        print("\n******************************************")

    make_change = force or yes_or_no("Make changes to {}?".format(filename))
    if make_change:
        with open(filename, "w") as f:
            f.write(modified)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="stella input file auto-upgrader")

    parser.add_argument("files", action="store", nargs="+", help="Input files")

    force_patch_group = parser.add_mutually_exclusive_group()
    force_patch_group.add_argument(
        "--force", "-f", action="store_true", help="Make changes without asking"
    )
    force_patch_group.add_argument(
        "--patch-only", "-p", action="store_true", help="Print the patches and exit"
    )

    parser.add_argument(
        "--quiet", "-q", action="store_true", help="Don't print patches"
    )

    args = parser.parse_args()

    for filename in args.files:
        modified = apply_fixes(filename, REPLACEMENTS)

        if modified is None:
            if not args.quiet:
                print("No changes to make to {}".format(filename))
            continue

        patch = create_patch(filename, modified)

        if args.patch_only:
            print(patch)
            continue

        possibly_apply_patch(filename, patch, modified, args.quiet, args.force)
