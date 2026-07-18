import type { MitoAnalysisData } from '@mito-architect/visualization-lib';

export const demoData: MitoAnalysisData = {
  metadata: {
    schema_version: '0.5',
    engine_version: '0.5.0-dev',
    sample: 'demo-mixed-haplogroups',
    input_path: 'fixtures/tiny.fastq',
    reference_path: 'core/data/rcrs.fasta',
    reference_accession: 'NC_012920.1',
    reference_length: 16569,
    threads: 1,
    resources: [
      {
        name: 'rCRS',
        version: 'NC_012920.1',
        path: 'rcrs.fasta',
        sha256: 'fc392cde8e63b4d2e3a870bb97cc0626dea33d46dfb8abdebffada040f42ec92',
        source: 'https://www.ncbi.nlm.nih.gov/nuccore/NC_012920.1',
        license: 'NCBI data terms',
        retrieved: 'bundled'
      },
      {
        name: 'phylotree-rcrs',
        version: '17.3',
        path: 'phylotree-rcrs-17.3.xml',
        sha256: '5490d3756503654001d07ed3a733bcfa8ee9ee225c50cd07c8c4260d239d85b1',
        source: 'https://github.com/genepi/phylotree-rcrs-17/releases/tag/17.3',
        license: 'MIT',
        retrieved: '2026-07-12'
      }
    ],
    algorithm_notes: [
      'FASTQ/SAM vertical slice',
      'NUMT specificity requires competitive nuclear-plus-mitochondrial alignment',
      'Cluster IDs are deterministic feature signatures until HDBSCAN adapter is enabled'
    ]
  },
  filter_stats: {
    input_reads: 7,
    passed_reads: 6,
    numt_filtered_reads: 1,
    numt_threshold: 0.3,
    input_alignment_records: 7,
    input_molecules: 7,
    numt_assessment: {
      mode: 'unaligned_fastq',
      nuclear_contigs_present: false,
      specificity_assessable: false
    }
  },
  coverage_metrics: {
    mean_depth: 13.333,
    pct_bins_gt20x: 5.556,
    max_depth: 22,
    mapping_quality_histogram: []
  },
  genes: [
    { name: 'MT-RNR1', start: 648, end: 1601, strand: '+', biotype: 'rRNA' },
    { name: 'MT-RNR2', start: 1671, end: 3229, strand: '+', biotype: 'rRNA' },
    { name: 'MT-ND1', start: 3307, end: 4262, strand: '+', biotype: 'protein_coding' },
    { name: 'MT-ND2', start: 4470, end: 5511, strand: '+', biotype: 'protein_coding' },
    { name: 'MT-CO1', start: 5904, end: 7445, strand: '+', biotype: 'protein_coding' },
    { name: 'MT-CO2', start: 7586, end: 8269, strand: '+', biotype: 'protein_coding' },
    { name: 'MT-ATP8', start: 8366, end: 8572, strand: '+', biotype: 'protein_coding' },
    { name: 'MT-ATP6', start: 8527, end: 9207, strand: '+', biotype: 'protein_coding' },
    { name: 'MT-CO3', start: 9207, end: 9990, strand: '+', biotype: 'protein_coding' },
    { name: 'MT-ND3', start: 10059, end: 10404, strand: '+', biotype: 'protein_coding' },
    { name: 'MT-ND4L', start: 10470, end: 10766, strand: '+', biotype: 'protein_coding' },
    { name: 'MT-ND4', start: 10760, end: 12137, strand: '+', biotype: 'protein_coding' },
    { name: 'MT-ND5', start: 12337, end: 14148, strand: '+', biotype: 'protein_coding' },
    { name: 'MT-ND6', start: 14149, end: 14673, strand: '-', biotype: 'protein_coding' },
    { name: 'MT-CYB', start: 14747, end: 15887, strand: '+', biotype: 'protein_coding' }
  ],
  coverage: Array.from({ length: 90 }, (_, index) => {
    const start = index * 184 + 1;
    return {
      start,
      end: Math.min(16569, start + 183),
      depth: 8 + Math.round(Math.sin(index / 6) * 4 + (index % 11))
    };
  }),
  svs: [
    {
      id: 'del:8483-13459',
      type: 'del',
      start: 8483,
      end: 13459,
      length: 4977,
      known_event: true,
      supporting_reads: ['demo-hapA-001 sv=del:8483-13459', 'demo-hapA-002 sv=del:8483-13459']
    }
  ],
  clusters: [
    {
      id: 0,
      label: 'H1',
      haplogroup: 'H',
      size: 3,
      consensus_haplotype: 'feature-consensus-0',
      reads: [
        'demo-hapA-001 sv=del:8483-13459',
        'demo-hapA-002 sv=del:8483-13459',
        'demo-hapA-003'
      ],
      sv_signature: [{ sv_id: 'del:8483-13459', support: 2 }]
    },
    {
      id: 1,
      label: 'H2',
      haplogroup: 'J',
      size: 2,
      consensus_haplotype: 'feature-consensus-1',
      reads: ['demo-hapB-001', 'demo-hapB-002'],
      sv_signature: []
    },
    {
      id: 2,
      label: 'H3',
      haplogroup: 'T',
      size: 1,
      consensus_haplotype: 'feature-consensus-2',
      reads: ['demo-hapC-001'],
      sv_signature: []
    }
  ],
  reads: [
    {
      id: 'demo-hapA-001 sv=del:8483-13459',
      length: 48,
      mean_quality: 40,
      numt_score: 0.05,
      filtered_numt: false,
      cluster_id: 0,
      snps: [
        { position: 73, ref: 'A', alt: 'G' },
        { position: 263, ref: 'A', alt: 'G' },
        { position: 8860, ref: 'A', alt: 'G' },
        {
          position: 8993,
          ref: 'T',
          alt: 'G',
          gene: 'MT-ATP6',
          consequence: 'missense_variant',
          protein: 'MT-ATP6',
          residue: 'p.Leu156Arg',
          annotation: {
            phenotype: 'NARP; Leigh syndrome',
            pathogenicity: 'pathogenic',
            references: ['PMID:1839049', 'PMID:1737708'],
            source: 'local-cache',
            sources: ['MITOMAP', 'ClinVar']
          },
          structure: {
            structure_id: '8H9S',
            chain: 'N',
            residue_index: 156,
            residue_label: 'p.Leu156Arg',
            complex: 'Complex V ATP synthase'
          }
        }
      ],
      sv_ids: ['del:8483-13459']
    },
    {
      id: 'demo-hapA-002 sv=del:8483-13459',
      length: 48,
      mean_quality: 39,
      numt_score: 0.05,
      filtered_numt: false,
      cluster_id: 0,
      snps: [
        { position: 73, ref: 'A', alt: 'G' },
        { position: 8860, ref: 'A', alt: 'G' },
        {
          position: 8993,
          ref: 'T',
          alt: 'G',
          gene: 'MT-ATP6',
          consequence: 'missense_variant',
          protein: 'MT-ATP6',
          residue: 'p.Leu156Arg',
          annotation: {
            phenotype: 'NARP; Leigh syndrome',
            pathogenicity: 'pathogenic',
            references: ['PMID:1839049', 'PMID:1737708'],
            source: 'local-cache',
            sources: ['MITOMAP', 'ClinVar']
          },
          structure: {
            structure_id: '8H9S',
            chain: 'N',
            residue_index: 156,
            residue_label: 'p.Leu156Arg',
            complex: 'Complex V ATP synthase'
          }
        },
        { position: 15326, ref: 'A', alt: 'G' }
      ],
      sv_ids: ['del:8483-13459']
    },
    {
      id: 'demo-hapA-003',
      length: 48,
      mean_quality: 38,
      numt_score: 0.05,
      filtered_numt: false,
      cluster_id: 0,
      snps: [
        { position: 73, ref: 'A', alt: 'G' },
        { position: 263, ref: 'A', alt: 'G' }
      ],
      sv_ids: []
    },
    {
      id: 'demo-hapB-001',
      length: 48,
      mean_quality: 40,
      numt_score: 0.05,
      filtered_numt: false,
      cluster_id: 1,
      snps: [
        { position: 750, ref: 'A', alt: 'G' },
        { position: 1438, ref: 'A', alt: 'G' },
        {
          position: 11778,
          ref: 'G',
          alt: 'A',
          gene: 'MT-ND4',
          consequence: 'missense_variant',
          protein: 'MT-ND4',
          residue: 'p.Arg340His',
          annotation: {
            phenotype: 'Leber hereditary optic neuropathy',
            pathogenicity: 'pathogenic',
            references: ['PMID:3413054', 'PMID:1679310'],
            source: 'local-cache',
            sources: ['MITOMAP', 'ClinVar']
          },
          structure: {
            structure_id: '9I4I',
            chain: 'r',
            residue_index: 340,
            residue_label: 'p.Arg340His',
            complex: 'Complex I NADH dehydrogenase'
          }
        }
      ],
      sv_ids: []
    },
    {
      id: 'demo-hapB-002',
      length: 48,
      mean_quality: 37,
      numt_score: 0.05,
      filtered_numt: false,
      cluster_id: 1,
      snps: [
        { position: 750, ref: 'A', alt: 'G' },
        {
          position: 11778,
          ref: 'G',
          alt: 'A',
          gene: 'MT-ND4',
          consequence: 'missense_variant',
          protein: 'MT-ND4',
          residue: 'p.Arg340His',
          annotation: {
            phenotype: 'Leber hereditary optic neuropathy',
            pathogenicity: 'pathogenic',
            references: ['PMID:3413054', 'PMID:1679310'],
            source: 'local-cache',
            sources: ['MITOMAP', 'ClinVar']
          },
          structure: {
            structure_id: '9I4I',
            chain: 'r',
            residue_index: 340,
            residue_label: 'p.Arg340His',
            complex: 'Complex I NADH dehydrogenase'
          }
        },
        { position: 3010, ref: 'G', alt: 'A' }
      ],
      sv_ids: []
    },
    {
      id: 'demo-hapC-001',
      length: 48,
      mean_quality: 40,
      numt_score: 0.05,
      filtered_numt: false,
      cluster_id: 2,
      snps: [
        { position: 2706, ref: 'A', alt: 'G' },
        {
          position: 3243,
          ref: 'A',
          alt: 'G',
          gene: 'MT-TL1',
          consequence: 'tRNA_variant',
          annotation: {
            phenotype: 'MELAS; MIDD',
            pathogenicity: 'pathogenic',
            references: ['PMID:2106746', 'PMID:8358449'],
            source: 'local-cache',
            sources: ['MITOMAP', 'ClinVar']
          }
        },
        { position: 7028, ref: 'C', alt: 'T' }
      ],
      sv_ids: []
    },
    {
      id: 'demo-numt nuclear',
      length: 48,
      mean_quality: 40,
      numt_score: 0.82,
      filtered_numt: true,
      cluster_id: -1,
      snps: [],
      sv_ids: []
    }
  ]
};
