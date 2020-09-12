# MAb Table Data Editing Instruction

## List of all fields

| Field         | Type    | Description                  | Required? |
| :---          | :---    | :---                         | :---      |
| `references`  | list    | List of `reference` objects. | Required  |
| `source`      | list    | List of free text source.    | Required  |
| `antibodies`  | list    | List of `antibody` objects.  | Required  |
| `config`      | text    | Configuration (free text).   | Required  |


### Fields of `reference` object

| Field                    | Type      | Description                                             | Required? |
| :---                     | :---      | :---                                                    | :---      |
| `doi`                    | text      | DOI number of a reference.                              | Required  |
| `firstAuthor`            | object    | First author.                                           | Optional  |
| `firstAuthor.surname`    | text      | First author's surname.                                 | Optional  |
| `firstAuthor.givenNames` | text      | First author's given names.                             | Optional  |
| `year`                   | number    | Publication year.                                       | Optional  |
| `journal`                | text      | Journal title.                                          | Optional  |
| `journalShort`           | text      | Journal short title (if not the same as journal title). | Optional  |

