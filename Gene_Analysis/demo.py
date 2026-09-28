"""Small, entirely fictional fixtures for an offline installation check."""
from pathlib import Path

GENOME = """# SYNTHETIC DEMO: not a person or real clinical evidence.
# genome build GRCh37
# strand plus
rsid\tchromosome\tposition\tgenotype
rs900000001\t1\t10001\tAG
rs900000002\t2\t10002\tCC
rs900000003\t3\t10003\tCT
"""
VCF = """##fileformat=VCFv4.2
##fileDate=20260901
##reference=GRCh37
##source=SYNTHETIC_DEMO_NOT_CLINICAL_EVIDENCE
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
1\t10001\t900000001\tA\tG\t.\t.\tGENEINFO=DEMO1:900000001;CLNSIG=Benign;CLNREVSTAT=reviewed_by_expert_panel;CLNDN=Fictional_demo_condition
2\t10002\t900000002\tC\tT\t.\t.\tGENEINFO=DEMO2:900000002;CLNSIG=Uncertain_significance;CLNREVSTAT=criteria_provided,_single_submitter;CLNDN=Fictional_demo_condition
3\t10003\t900000003\tC\tT\t.\t.\tGENEINFO=DEMO3:900000003;CLNSIG=Uncertain_significance;CLNREVSTAT=criteria_provided,_single_submitter;CLNDN=Fictional_demo_condition
"""
CLINGEN = '''"CLINGEN GENE DISEASE VALIDITY CURATIONS","","","","","","","","",""
"FILE CREATED: 2026-09-01","","","","","","","","",""
"WEBPAGE: SYNTHETIC DEMO","","","","","","","","",""
"+++++++++++","+++++++++++","+++++++++++","+++++++++++","+++++++++++","+++++++++++","+++++++++++","+++++++++++","+++++++++++","+++++++++++"
"GENE SYMBOL","GENE ID (HGNC)","DISEASE LABEL","DISEASE ID (MONDO)","MOI","SOP","CLASSIFICATION","ONLINE REPORT","CLASSIFICATION DATE","GCEP"
"DEMO1","HGNC:900000001","Fictional demo condition","MONDO:9000001","AD","SOP10","Limited","https://example.invalid/demo","2026-09-01","Fictional demo panel"
'''


def run_demo(output_dir: str) -> int:
    from .cli import build_parser, cmd_analyze
    from .evidence.clinical import refresh_clinvar_database
    from .evidence.gene_validity import refresh_gene_validity
    output = Path(output_dir).expanduser().resolve()
    # Keep the tiny pinned source/database files with the report for replay.
    fixtures = output / 'synthetic_evidence'
    fixtures.mkdir(parents=True, exist_ok=True)
    genome = fixtures / 'synthetic.txt'
    vcf = fixtures / 'synthetic.vcf'
    genes = fixtures / 'synthetic_clingen.csv'
    for path, text in [(genome, GENOME), (vcf, VCF), (genes, CLINGEN)]:
        path.write_text(text, encoding='utf-8')
    manifest = refresh_clinvar_database(cache_dir=fixtures, source_path=vcf)
    refresh_gene_validity(cache_dir=fixtures, source_path=genes)
    args = build_parser().parse_args(['analyze', str(genome), '--offline',
        '--clinvar-manifest', manifest['manifest_path'], '--cache-dir', str(fixtures),
        '--output-dir', str(output)])
    args.synthetic_demo = True
    print('SYNTHETIC DEMO — fictional variants and classifications; no medical meaning.')
    return cmd_analyze(args)
