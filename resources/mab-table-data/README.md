# MAb Table Data Editing Instruction

Please refer to `mab-table-data/` folder for examples.

## List of all fields

| Field         | Type    | Description                  | Required? |
| :---          | :---    | :---                         | :---      |
| `references`  | list    | List of `Reference` objects. | Required  |
| `source`      | list    | List of free text source.    | Required  |
| `antibodies`  | list    | List of `Antibody` objects.  | Required  |
| `config`      | text    | Configuration (free text).   | Required  |


### `Reference` object

| Field                    | Type      | Description                                             | Required? |
| :---                     | :---      | :---                                                    | :---      |
| `doi`                    | text      | DOI number of a reference.                              | Required  |
| `firstAuthor`            | object    | First author.                                           | Optional  |
| `firstAuthor.surname`    | text      | First author's surname.                                 | Optional  |
| `firstAuthor.givenNames` | text      | First author's given names.                             | Optional  |
| `year`                   | number    | Publication year.                                       | Optional  |
| `journal`                | text      | Journal title.                                          | Optional  |
| `journalShort`           | text      | Journal short title (if not the same as journal title). | Optional  |

### `Antibody` object

| Field                       | Type        | Description                                                                  | Required? |
| :---                        | :---        | :---                                                                         | :---      |
| name                        | text        | Name of this MAb.                                                            | Required  |
| ec50                        | text/number | The EC50 info of this MAb. If only puts in number the default unit is ng/ml. | Optional  |
| ec50Note                    | text        | Notation for EC50. Can be "PV", "IC100" or "IC90".                           | Optional  |
| pdb                         | text        | Descriptional text for PDB structures.                                       | Optional  |
| species                     | text        |                                                                              | Optional  |
| IGHV/PcntMutH/CDRH3Len/IGHJ | -           | Auto populated fields from MAb alignments. Don't change them.                | Optional  |
| IGLV/PcntMutL/CDRL3Len/IGLJ | -           | Auto populated fields from MAb alignments. Don't change them.                | Optional  |
