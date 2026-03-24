
import pandas as pd

df = pd.read_csv('/mnt/c/home/git/simtiva/data/supplementary_digital_content_3.txt', sep=r'\s+')

lines = []
lines.append(f"Shape: {df.shape}")
lines.append(f"Columns: {list(df.columns)}")
lines.append("")
lines.append(f"DVTY unique: {sorted(df['DVTY'].unique().tolist())}")
lines.append(f"M1F2 unique: {sorted(df['M1F2'].unique().tolist())}")
lines.append(f"A1V2 unique: {sorted(df['A1V2'].unique().tolist())}")
lines.append(f"P1V2 unique: {sorted(df['P1V2'].unique().tolist())}")
lines.append(f"STDY unique: {sorted(df['STDY'].unique().tolist())}")
lines.append(f"GRP unique: {sorted(df['GRP'].unique().tolist())}")
lines.append(f"TECH unique: {sorted(df['TECH'].unique().tolist())}")
lines.append(f"EVID unique: {sorted(df['EVID'].unique().tolist())}")
lines.append("")
lines.append(f"Unique ID (subjects): {df['ID'].nunique()}")
lines.append(f"Unique DID: {df['DID'].nunique()}")
lines.append("")

obs = df[df['EVID'] == 0]
dose = df[df['EVID'] == 1]
lines.append(f"Total rows: {len(df)}")
lines.append(f"Dosing events (EVID=1): {len(dose)}")
lines.append(f"Observations (EVID=0): {len(obs)}")
lines.append(f"  DVTY breakdown in obs: {obs['DVTY'].value_counts().to_dict()}")
obs_nonzero = obs[obs['DV'] != 0]
lines.append(f"  Non-zero DV rows: {len(obs_nonzero)}")
neg_dv = obs[obs['DV'] < 0]
lines.append(f"  Negative DV rows (LLOQ/BLQ?): {len(neg_dv)}")
lines.append("")
lines.append(f"Age range: {df['AGE'].min()} - {df['AGE'].max()}")
lines.append(f"Weight range: {df['WGT'].min()} - {df['WGT'].max()}")
lines.append(f"Height range (HGT): {df['HGT'].min()} - {df['HGT'].max()}")
lines.append(f"PMA range: {df['PMA'].min()} - {df['PMA'].max()}")
lines.append(f"FFM range (FFMZ): {df['FFMZ'].min()} - {df['FFMZ'].max()}")
lines.append(f"BMI range: {df['BMI'].min()} - {df['BMI'].max()}")
lines.append("")
lines.append(f"AMT stats (dosing): ")
lines.append(str(df[df['EVID']==1]['AMT'].describe()))
lines.append("")
lines.append(f"RATE stats (dosing): ")
lines.append(str(df[df['EVID']==1]['RATE'].describe()))
lines.append("")
lines.append(f"DV stats (non-zero obs): ")
lines.append(str(obs_nonzero['DV'].describe()))
lines.append("")
lines.append("Sample of first few non-zero DV rows:")
lines.append(obs_nonzero[['ID','TIME','DV','DVTY','AGE','WGT','M1F2']].head(10).to_string())

print('\n'.join(lines))
