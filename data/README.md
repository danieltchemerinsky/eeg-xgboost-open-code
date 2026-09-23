# Data

Place the two data files here before running the analysis:

```
data/outputminus.RData
data/AllPatientsEEGandETCAbsentObs5.RData
```

Both are saved R workspaces. The scripts load them into a private environment
and take only the object `combined`, because these files contain other objects
that would otherwise overwrite analysis variables.

| File | Contents | Used by |
|---|---|---|
| `outputminus.RData` | 121 patients × 20 columns: `Outcome`, `Patient`, and 18 EEG reactivity features (delta/theta/alpha × six stimulation sites) | `01`, `09` |
| `AllPatientsEEGandETCAbsentObs5.RData` | The same cohort with four additional clinical variables: `shockable_rhythm`, `basic_cpr`, `SEP`, `age` | `02`, `05` |

Patients with more than five missing observations are excluded, leaving
N = 114 (CPC 1: 59, CPC 2: 22, CPC 3: 2, CPC 5: 31; no CPC 4 observed).

**These files are not distributed with this repository, and will not be.** They
come from the TTH48 multicentre trial (Kirkegaard et al., 2017) and contain
clinical data on comatose cardiac arrest patients. Data-sharing agreements
between the participating centres and European data-protection law, including
the GDPR, do not permit open publication, so they cannot be deposited in this or
any other public repository.

De-identified data may be requested from the corresponding author, subject to
approval by the TTH48 Trial Steering Committee and a data-access agreement under
the host institutions' governance framework. See the data availability statement
in the paper.

This directory is empty by design apart from this file. `.gitignore` excludes
`*.RData` and `output/`, so the repository holds only code.
