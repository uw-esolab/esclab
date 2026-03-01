# esclab
Energy systems simulation and control

## Conversion queue automation

Generate a prioritized Fortran-to-Python conversion queue with:

`python scripts/generate_conversion_queue.py`

This writes:

- `docs/conversion_queue.json` (machine-readable queue)
- `docs/conversion_queue.md` (human-readable backlog)

Optional:

- `python scripts/generate_conversion_queue.py --print-top 25`
