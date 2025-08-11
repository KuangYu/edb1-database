## Analysis of EDB1 Database

### Key scripts:
`analysis.py`: prescreening to find out target mols/salts, and count the coverage of formulas by target mols/salts
`typification.py`: analyze the atomtypes in target mols/salts, find out representitive ones for experiments
`query_purchase.py`: query how fast the chemicals in `mols/salts_selected.dat` can be purchased?

### Key files:
`mols_target.dat`: target mols/salts
`atomtypes.dat`: involved atomtypes
`mols_not_covered.dat`: mols with atomtypes not covered by selected molecules, requires manual handling
`mols_selected.dat`: mols selected for experiments
`salts_selected.dat`: salts selected for experiments
`index_formula.txt`: covered formula indices

### External (input) csv data files:
`edb1.csv`: EDB1 formula library
`prices.csv`: Availability and prices of chemicals
`order_info.csv`: Order availability

