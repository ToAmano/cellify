# Example 5: Retrieve Crystal Structures by Chemical Formula

This example demonstrates how to use `cellify` to retrieve crystal structures from public databases (Materials Project and Crystallography Open Database) by providing a chemical formula.

## Description

When the input parameter `-i`/`--input` is a chemical formula (instead of a path to an existing file), `cellify` queries public OPTIMADE servers to search for corresponding crystal structures.

- **Si**: Shows how `cellify` automatically parses the formula and queries both databases, returning polymorphic cell sizes like standard cubic diamond silicon.
- **TiO2**: Illustrates that `cellify` handles composition order automatically (Hill notation `O2Ti`), returning all polymorphs of titanium dioxide (e.g., rutile, anatase) across different unit cell representations.
- **H3S**: Demonstrates querying a computationally-discovered high-pressure phase which is present in DFT databases (Materials Project) but not in experimental repositories (COD).

## How to Run

Execute the provided script to run the queries:

```bash
bash run.sh
```
