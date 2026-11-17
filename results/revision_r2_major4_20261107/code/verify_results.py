from pathlib import Path
import json,hashlib
import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_rand_score,normalized_mutual_info_score
ROOT=Path('/home/xd/y/Benchmark');out=ROOT/'results/revision_r2_major4_20261107'
rows=pd.read_csv(out/'section_seed_metrics.tsv',sep='\t',dtype={'heldout_slice':str})
assert len(rows)==72 and not rows.duplicated(['heldout_donor','heldout_slice','strategy','condition','seed']).any()
counts=pd.read_csv(out/'selection_audit.tsv',sep='\t')
mapping=pd.read_csv(out/'target_mapping.tsv',sep='\t').fillna('')
for r in counts.itertuples():
    genes=json.loads((out/'folds'/r.heldout_donor/(r.strategy+('_retain' if r.condition=='Retain targets' else '_base')+'_genes.json')).read_text())
    targets=set(mapping.loc[mapping.heldout_donor.eq(r.heldout_donor),'mapped_gene'])-{''}
    assert len(set(genes))==len(genes)==r.actual_features==r.requested_features==2000
    assert len(targets)==r.available_targets and len(targets&set(genes))==r.retained_targets
    if r.condition=='Retain targets':assert targets.issubset(genes)
for r in rows.itertuples():
    tag=r.strategy+('_retain' if r.condition=='Retain targets' else '_base')
    p=np.load(out/'folds'/r.heldout_donor/'predictions'/f'{tag}__{r.heldout_slice}__seed{r.seed}.npz')
    assert abs(adjusted_rand_score(p['truth'][p['valid']],p['pred'][p['valid']])-r.heldout_ari)<1e-12
    assert abs(normalized_mutual_info_score(p['truth'][p['valid']],p['pred'][p['valid']])-r.heldout_nmi)<1e-12
hashes=json.loads((out/'input_hashes.json').read_text())
assert all(hashlib.sha256(Path(p).read_bytes()).hexdigest()==h for p,h in hashes.items())
audit={'all_72_predictions_independently_rescored':True,'all_12_feature_sets_checked':True,'exact_budget_2000':True,'forced_sets_contain_all_available_targets':True,'prior_inputs_unchanged':True}
(out/'independent_verification.json').write_text(json.dumps(audit,indent=2));print(json.dumps(audit))
print(pd.read_csv(out/'mean_summary.tsv',sep='\t').to_string(index=False))
print(pd.read_csv(out/'paired_donor_differences.tsv',sep='\t')[['heldout_donor','strategy','delta_heldout_ari','delta_heldout_nmi','delta_heldout_moran_percentile']].to_string(index=False))
