import type {
  AggregateVariant,
  AlignmentFragment,
  ClinicalAnnotation,
  EvidenceObservation,
  MitoAnalysisData,
  ProteinStructureMapping,
  ReadFeature
} from '@mito-architect/visualization-lib';
import { useQuery } from '@tanstack/react-query';
import {
  Activity,
  Atom,
  ChartNoAxesCombined,
  Database,
  Dna,
  GitFork,
} from 'lucide-react';
import { type KeyboardEvent, useEffect, useMemo, useRef, useState } from 'react';
import { FixedSizeList, type ListChildComponentProps } from 'react-window';
import {
  type EvidenceSearchFilters,
  type EvidenceSearchRow,
  getEvidencePage,
  searchEvidence
} from '../lib/api';
import { nextTabIndex, normalizeExactEvidenceFilter } from '../lib/analysisNavigation';
import { useMitoStore } from '../lib/store';
import NglProteinViewer, { type ProteinViewerVariant } from './NglProteinViewer';

interface AnalysisModulesPanelProps {
  data: MitoAnalysisData;
  jobId?: string;
}

type ModuleTab =
  | 'haplogroups'
  | 'clinical'
  | 'coverage'
  | 'evidence'
  | 'rearrangements'
  | 'protein';

interface VariantSummary {
  key: string;
  label: string;
  type: string;
  position: number;
  support: number;
  frequency: number;
  gene?: string;
  consequence?: string;
  protein?: string;
  residue?: string;
  annotation?: ClinicalAnnotation;
  structure?: ProteinStructureMapping;
  molecules: Array<{ id: string; clusterId: number }>;
  callableDepth?: number;
  ci95Low?: number;
  ci95High?: number;
  strandSupport?: AggregateVariant['strand_support'];
  strandBiasDelta?: number | null;
  alleleQuality?: AggregateVariant['allele_quality'];
  readPosition?: AggregateVariant['read_position'];
}

const MODULES: Array<{ id: ModuleTab; label: string; Icon: typeof Dna }> = [
  { id: 'haplogroups', label: 'Haplogroups', Icon: Dna },
  { id: 'clinical', label: 'Clinical', Icon: Database },
  { id: 'coverage', label: 'Coverage', Icon: Activity },
  { id: 'evidence', label: 'Molecules & phase', Icon: GitFork },
  { id: 'rearrangements', label: 'Rearrangements', Icon: GitFork },
  { id: 'protein', label: '3D Protein', Icon: Atom }
];

export default function AnalysisModulesPanel({ data, jobId }: AnalysisModulesPanelProps) {
  const [activeTab, setActiveTab] = useState<ModuleTab>('haplogroups');
  const tabRefs = useRef<Array<HTMLButtonElement | null>>([]);
  const passedReads = useMemo(() => data.reads.filter((read) => !read.filtered_numt), [data.reads]);
  const variants = useMemo(() => summarizeVariants(passedReads, data), [passedReads, data]);
  const coverageMetrics = useMemo(() => summarizeCoverage(data), [data]);
  const active = MODULES.find((module) => module.id === activeTab) ?? MODULES[0];
  const handleTabKeyDown = (event: KeyboardEvent<HTMLButtonElement>, index: number) => {
    const nextIndex = nextTabIndex(event.key, index, MODULES.length);
    if (nextIndex === undefined) return;
    event.preventDefault();
    setActiveTab(MODULES[nextIndex].id);
    tabRefs.current[nextIndex]?.focus();
  };

  return (
    <section className="rounded-lg border border-line bg-panel">
      <div className="flex flex-wrap items-center justify-between gap-3 border-b border-line px-4 py-3">
        <h2 className="flex items-center gap-2 text-sm font-semibold tracking-normal">
          <ChartNoAxesCombined className="h-4 w-4 text-aqua" aria-hidden />
          Analysis Modules
        </h2>
        <div className="flex flex-wrap gap-2" role="tablist" aria-label="Analysis modules">
          {MODULES.map(({ id, label, Icon }, index) => (
            <button
              key={id}
              id={`analysis-tab-${id}`}
              type="button"
              role="tab"
              aria-selected={activeTab === id}
              aria-controls={`analysis-panel-${id}`}
              tabIndex={activeTab === id ? 0 : -1}
              ref={(element) => { tabRefs.current[index] = element; }}
              onClick={() => setActiveTab(id)}
              onKeyDown={(event) => handleTabKeyDown(event, index)}
              className={[
                'inline-flex items-center gap-2 rounded-md border px-3 py-2 text-sm',
                activeTab === id ? 'border-aqua text-aqua' : 'border-line text-muted hover:bg-panel2'
              ].join(' ')}
            >
              <Icon className="h-4 w-4" aria-hidden />
              {label}
            </button>
          ))}
        </div>
      </div>

      <div
        id={`analysis-panel-${activeTab}`}
        role="tabpanel"
        aria-labelledby={`analysis-tab-${activeTab}`}
        className="p-4"
      >
        <div className="mb-4 flex items-center gap-2 text-sm font-semibold">
          <active.Icon className="h-4 w-4 text-aqua" aria-hidden />
          {active.label}
        </div>
        {activeTab === 'haplogroups' && <HaplogroupPanel data={data} />}
        {activeTab === 'clinical' && <ClinicalPanel variants={variants} />}
        {activeTab === 'coverage' && <CoveragePanel data={data} metrics={coverageMetrics} />}
        {activeTab === 'evidence' && <EvidencePanel data={data} jobId={jobId} />}
        {activeTab === 'rearrangements' && <RearrangementsPanel data={data} />}
        {activeTab === 'protein' && <ProteinPanel variants={variants} />}
      </div>
    </section>
  );
}

function EvidencePanel({ data, jobId }: { data: MitoAnalysisData; jobId?: string }) {
  const callability = data.callability ?? [];
  const events = data.events ?? [];
  const phaseLinks = data.phase_links ?? [];
  const remoteMoleculeIds = useMemo(
    () => (data.molecules ?? []).slice(0, 500).map((molecule) => molecule.id),
    [data.molecules]
  );
  const remoteEventIds = useMemo(
    () => (data.events ?? []).slice(0, 500).map((event) => event.id),
    [data.events]
  );
  const moleculesById = new Map((data.molecules ?? []).map((molecule) => [molecule.id, molecule]));
  const selectedEventId = useMitoStore((state) => state.selectedEventId);
  const selectedMoleculeId = useMitoStore((state) => state.selectedMoleculeId);
  const selectedPhaseId = useMitoStore((state) => state.selectedPhaseId);
  const setSelectedEventId = useMitoStore((state) => state.setSelectedEventId);
  const setSelectedMoleculeId = useMitoStore((state) => state.setSelectedMoleculeId);
  const setSelectedPhaseId = useMitoStore((state) => state.setSelectedPhaseId);
  const remotePages = data.evidence_encoding?.observation_storage === 'remote_http_pages';
  const pageCount = data.evidence_encoding?.observation_page_count ?? 0;
  const [remotePageIndex, setRemotePageIndex] = useState(0);
  const remotePage = useQuery({
    queryKey: ['evidence-page', jobId, remotePageIndex],
    queryFn: ({ signal }) => getEvidencePage(jobId!, remotePageIndex, signal),
    enabled: remotePages && Boolean(jobId) && pageCount > 0,
    staleTime: Number.POSITIVE_INFINITY,
    gcTime: 5 * 60 * 1000
  });
  useEffect(() => {
    if (remotePageIndex >= pageCount) setRemotePageIndex(Math.max(0, pageCount - 1));
  }, [pageCount, remotePageIndex]);
  const selection = selectedPhaseId
    ? { kind: 'phase' as const, id: selectedPhaseId }
    : selectedEventId
      ? { kind: 'event' as const, id: selectedEventId }
      : selectedMoleculeId
        ? { kind: 'molecule' as const, id: selectedMoleculeId }
        : null;
  const selected = selection && (
    selection.kind === 'event'
      ? events.some((event) => event.id === selection.id)
      : selection.kind === 'molecule'
        ? callability.some((item) => item.molecule_id === selection.id)
        : phaseLinks.some((link) => link.id === selection.id)
  ) ? selection : null;
  const selectedRemoteObservations = useQuery({
    queryKey: ['evidence-selection', jobId, selected?.kind, selected?.id],
    queryFn: async ({ signal }) => {
      if (!selected || !jobId) return [];
      if (selected.kind === 'event') {
        return (await searchEvidence(jobId, { eventId: selected.id }, 0, 300, signal)).rows;
      }
      if (selected.kind === 'molecule') {
        return (await searchEvidence(jobId, { moleculeId: selected.id }, 0, 300, signal)).rows;
      }
      const phase = phaseLinks.find((link) => link.id === selected.id);
      if (!phase) return [];
      const responses = await Promise.all([
        searchEvidence(jobId, { eventId: phase.event_a_id }, 0, 300, signal),
        searchEvidence(jobId, { eventId: phase.event_b_id }, 0, 300, signal)
      ]);
      return Array.from(
        new Map(responses.flatMap((response) => response.rows).map((row) => [row.id, row])).values()
      )
        .sort((left, right) => left.page_index - right.page_index || left.row_index - right.row_index)
        .slice(0, 300);
    },
    enabled: remotePages && Boolean(jobId) && selected !== null,
    staleTime: Number.POSITIVE_INFINITY,
    gcTime: 5 * 60 * 1000
  });
  const selectedObservations = selected
    ? remotePages
      ? selectedRemoteObservations.data ?? []
      : matchingObservations(data, selected, 300)
    : [];
  const observationCount = data.evidence_encoding?.observation_count ?? data.observations?.length ?? 0;
  if (data.metadata.schema_version !== '0.6') {
    return (
      <div className="rounded-md border border-amber/50 bg-amber/10 p-4 text-sm text-amber">
        Molecule evidence is not present in this schema {data.metadata.schema_version ?? 'unknown'} result.
        Run analysis with the opt-in schema 0.6 evidence graph; a missing projection is not evidence that
        no molecular linkage exists.
      </div>
    );
  }

  return (
    <div className="grid gap-4">
      <div className="grid gap-3 sm:grid-cols-2 xl:grid-cols-4">
        <Metric label="Assembled molecules" value={(data.molecules?.length ?? 0).toString()} />
        <Metric label="Normalized events" value={events.length.toString()} />
        <Metric label="Phase links" value={phaseLinks.length.toString()} />
        <Metric label="Evidence rows" value={observationCount.toString()} />
      </div>

      {remotePages && (
        <>
          {jobId ? (
            <RemoteEvidenceSearch
              jobId={jobId}
              moleculeIds={remoteMoleculeIds}
              eventIds={remoteEventIds}
            />
          ) : (
            <div role="alert" className="rounded border border-amber/50 bg-amber/10 px-3 py-2 text-sm text-amber">
              Remote evidence search requires the originating job identifier.
            </div>
          )}
          <RemoteEvidenceBrowser
            page={remotePage.data}
            pageIndex={remotePageIndex}
            pageCount={pageCount}
            loading={remotePage.isLoading}
            error={remotePage.error as Error | null}
            onPageChange={setRemotePageIndex}
          />
        </>
      )}

      <div className="overflow-hidden rounded-md border border-line">
        <div className="grid grid-cols-[1.05fr_0.7fr_0.65fr_1.1fr_1.4fr] gap-2 border-b border-line bg-panel2 px-3 py-2 text-xs font-semibold uppercase tracking-normal text-muted">
          <div>Molecule</div><div>Status</div><div>Callable</div><div>Protocol decision</div><div>Alignment provenance</div>
        </div>
        {callability.length === 0 ? <EmptyRow text="Callability projection is empty" /> : callability.slice(0, 100).map((item) => {
          const molecule = moleculesById.get(item.molecule_id);
          const protocolDecision = [
            molecule?.identity_policy,
            ...(molecule?.protocol_flags ?? []),
            ...(molecule?.exclusion_reasons ?? []).map((reason) => `excluded:${reason}`)
          ].filter(Boolean).join(', ');
          return (
            <button
              key={item.molecule_id}
              type="button"
              aria-pressed={selected?.kind === 'molecule' && selected.id === item.molecule_id}
              onClick={() => {
                setSelectedEventId(undefined);
                setSelectedPhaseId(undefined);
                setSelectedMoleculeId(item.molecule_id);
              }}
              className="grid w-full grid-cols-[1.05fr_0.7fr_0.65fr_1.1fr_1.4fr] gap-2 border-b border-line px-3 py-2 text-left text-sm last:border-0 hover:bg-panel2 aria-pressed:bg-panel2"
            >
              <div className="truncate font-semibold" title={item.molecule_id}>{item.molecule_id}</div>
              <div className={molecule?.evidence_eligible && item.known ? 'text-aqua' : 'text-amber'}>{item.status}</div>
              <div>{item.known ? formatPercent(item.callable_fraction) : 'not assessable'}</div>
              <div className="truncate text-muted" title={protocolDecision || 'no protocol metadata'}>
                {protocolDecision || 'default identity'}
              </div>
              <div className="truncate text-muted" title={item.alignments.map((alignment) => `${alignment.alignment_id}:${alignment.status}`).join(', ')}>
                {item.alignments.map((alignment) => `${alignment.alignment_id} ${alignment.status}`).join(', ') || 'no eligible alignment'}
              </div>
            </button>
          );
        })}
      </div>

      <div className="overflow-hidden rounded-md border border-line">
        <div className="grid grid-cols-[1.35fr_0.8fr_0.65fr_0.65fr_0.95fr] gap-2 border-b border-line bg-panel2 px-3 py-2 text-xs font-semibold uppercase tracking-normal text-muted">
          <div>Event</div><div>Type</div><div>ALT</div><div>Callable</div><div>Assessability</div>
        </div>
        {events.length === 0 ? <EmptyRow text="No normalized events" /> : events.slice(0, 200).map((event) => (
          <button
            key={event.id}
            type="button"
            aria-pressed={selected?.kind === 'event' && selected.id === event.id}
            onClick={() => {
              setSelectedMoleculeId(undefined);
              setSelectedPhaseId(undefined);
              setSelectedEventId(event.id);
            }}
            className="grid w-full grid-cols-[1.35fr_0.8fr_0.65fr_0.65fr_0.95fr] gap-2 border-b border-line px-3 py-2 text-left text-sm last:border-0 hover:bg-panel2 aria-pressed:bg-panel2"
          >
            <div className="truncate font-semibold" title={event.id}>{eventLabel(event)}</div>
            <div>{event.type}</div>
            <div>{event.evidence_counts.alternate}</div>
            <div>{event.evidence_counts.callable}</div>
            <div
              className={event.assessability === 'REFERENCE_AND_ALTERNATE' ? 'text-aqua' : 'text-amber'}
              title={event.negative_evidence_rule}
            >
              {event.assessability === 'REFERENCE_AND_ALTERNATE' ? 'complete' : 'support-only'}
            </div>
          </button>
        ))}
      </div>

      <div className="overflow-hidden rounded-md border border-line">
        <div className="grid grid-cols-[1fr_1fr_0.48fr_0.48fr_0.48fr_0.48fr_0.65fr_0.8fr] gap-2 border-b border-line bg-panel2 px-3 py-2 text-xs font-semibold uppercase tracking-normal text-muted">
          <div>Event A</div><div>Event B</div><div>Joint</div><div>A+/B−</div><div>A−/B+</div><div>Both ALT</div><div>Co-ALT CI</div><div>Inspect</div>
        </div>
        {phaseLinks.length === 0 ? <EmptyRow text="No phase links with co-observed alternate evidence" /> : phaseLinks.slice(0, 200).map((link) => (
          <div key={link.id} className="grid grid-cols-[1fr_1fr_0.48fr_0.48fr_0.48fr_0.48fr_0.65fr_0.8fr] gap-2 border-b border-line px-3 py-2 text-sm last:border-0">
            <button type="button" className="truncate text-left hover:text-aqua" title={link.event_a_id} onClick={() => { setSelectedMoleculeId(undefined); setSelectedPhaseId(undefined); setSelectedEventId(link.event_a_id); }}>{link.event_a_id}</button>
            <button type="button" className="truncate text-left hover:text-aqua" title={link.event_b_id} onClick={() => { setSelectedMoleculeId(undefined); setSelectedPhaseId(undefined); setSelectedEventId(link.event_b_id); }}>{link.event_b_id}</button>
            <div>{link.jointly_callable}</div>
            <div>{link.a_alternate_b_absent}</div>
            <div>{link.a_absent_b_alternate}</div>
            <div title={phaseMoleculeTitle(link.supporting_molecule_indices, data)}>{link.both_alternate}</div>
            <div title={`uncertain n=${link.jointly_uncertain}; delta=${link.linkage_delta.toFixed(4)}`}>{formatPercent(link.co_alternate_ci95_low)}–{formatPercent(link.co_alternate_ci95_high)}</div>
            <button type="button" aria-pressed={selected?.kind === 'phase' && selected.id === link.id} className={link.assessability === 'COMPLETE_FOR_BOTH_EVENTS' ? 'text-left text-aqua hover:underline' : 'text-left text-amber hover:underline'} onClick={() => { setSelectedEventId(undefined); setSelectedMoleculeId(undefined); setSelectedPhaseId(link.id); }}>
              {link.assessability === 'COMPLETE_FOR_BOTH_EVENTS' ? 'complete pair' : 'conditioned'}
            </button>
          </div>
        ))}
      </div>
      <div className="overflow-hidden rounded-md border border-line">
        <div className="flex items-center justify-between gap-2 border-b border-line bg-panel2 px-3 py-2 text-xs font-semibold uppercase tracking-normal text-muted">
          <span>{selected ? `${selected.kind}: ${selected.id}` : 'Evidence inspector'}</span>
          {selected && <button type="button" className="text-aqua hover:underline" onClick={() => { setSelectedEventId(undefined); setSelectedMoleculeId(undefined); setSelectedPhaseId(undefined); }}>clear</button>}
        </div>
        {!selected ? <EmptyRow text="Select a molecule, event, or phase endpoint to trace its evidence" /> : selectedRemoteObservations.isLoading ? (
          <EmptyRow text="Searching all immutable evidence pages…" />
        ) : selectedRemoteObservations.isError ? (
          <div role="alert" className="px-3 py-4 text-sm text-red-300">
            {(selectedRemoteObservations.error as Error).message}
          </div>
        ) : selectedObservations.length === 0 ? (
          <EmptyRow text={remotePages ? 'No explicit observation in the global evidence index; sparse absence is not REF' : 'No explicit observation; sparse absence means NOT_CALLABLE'} />
        ) : (
          <>
            <div className="grid grid-cols-[1fr_1fr_0.65fr_0.8fr_1fr] gap-2 border-b border-line px-3 py-2 text-xs font-semibold uppercase tracking-normal text-muted">
              <div>Molecule</div><div>Event</div><div>State</div><div>Alignment</div><div>Evidence</div>
            </div>
            {selectedObservations.slice(0, 300).map((observation) => (
              <div key={observation.id} className="grid grid-cols-[1fr_1fr_0.65fr_0.8fr_1fr] gap-2 border-b border-line px-3 py-2 text-sm last:border-0">
                <button type="button" className="truncate text-left hover:text-aqua" title={observation.molecule_id} onClick={() => { setSelectedEventId(undefined); setSelectedPhaseId(undefined); setSelectedMoleculeId(observation.molecule_id); }}>{observation.molecule_id}</button>
                <button type="button" className="truncate text-left hover:text-aqua" title={observation.event_id} onClick={() => { setSelectedMoleculeId(undefined); setSelectedPhaseId(undefined); setSelectedEventId(observation.event_id); }}>{observation.event_id}</button>
                <div>{observation.state}</div>
                <div className="truncate" title={observation.alignment_id}>{observation.alignment_id}</div>
                <div className="truncate text-muted" title={observation.evidence_source}>{observation.evidence_source}</div>
              </div>
            ))}
          </>
        )}
      </div>
      {(callability.length > 100 || events.length > 200 || phaseLinks.length > 200 || observationCount > 300) && (
        <div className="text-xs text-muted">Preview is bounded in the browser; observations remain in deterministic columnar pages.</div>
      )}
    </div>
  );
}

function eventLabel(event: NonNullable<MitoAnalysisData['events']>[number]) {
  if (event.type === 'SNV' && event.start && event.ref && event.alt) {
    return `m.${event.start}${event.ref}>${event.alt}`;
  }
  if (event.start && event.end) {
    return `${event.start}-${event.end}`;
  }
  return event.id;
}

function phaseMoleculeTitle(indices: number[], data: MitoAnalysisData): string {
  if (indices.length === 0) return 'no co-ALT support';
  const molecules = data.molecules ?? [];
  return indices.map((index) => molecules[index]?.id ?? `unresolved molecule index ${index}`).join(', ');
}

function matchingObservations(
  data: MitoAnalysisData,
  selected: { kind: 'event' | 'molecule' | 'phase'; id: string },
  limit: number,
  pageOverride?: NonNullable<MitoAnalysisData['observation_pages']>
): EvidenceObservation[] {
  const matches: EvidenceObservation[] = [];
  const phase = selected.kind === 'phase'
    ? (data.phase_links ?? []).find((link) => link.id === selected.id)
    : undefined;
  const accept = (observation: EvidenceObservation) => {
    const matched = selected.kind === 'event'
      ? observation.event_id === selected.id
      : selected.kind === 'molecule'
        ? observation.molecule_id === selected.id
        : phase !== undefined && (observation.event_id === phase.event_a_id || observation.event_id === phase.event_b_id);
    if (matched && matches.length < limit) matches.push(observation);
  };
  if (data.observations) {
    for (const observation of data.observations) {
      accept(observation);
      if (matches.length === limit) return matches;
    }
    return matches;
  }
  for (const page of pageOverride ?? data.observation_pages ?? []) {
    for (let index = 0; index < page.count; index += 1) {
      accept({
        id: `observation:${page.offset + index}`,
        molecule_id: page.columns.molecule_id[index],
        event_id: page.columns.event_id[index],
        alignment_id: page.columns.alignment_id[index],
        state: page.columns.state[index],
        observed_allele: page.columns.observed_allele[index],
        base_quality: page.columns.base_quality[index],
        mapping_quality: page.columns.mapping_quality[index],
        strand: page.columns.strand[index],
        evidence_source: page.columns.evidence_source[index],
        read_position: page.columns.read_position[index]
      });
      if (matches.length === limit) return matches;
    }
  }
  return matches;
}

const REMOTE_SEARCH_LIMIT = 100;

function RemoteEvidenceSearch({
  jobId,
  moleculeIds,
  eventIds
}: {
  jobId: string;
  moleculeIds: string[];
  eventIds: string[];
}) {
  const moleculeSuggestions = useMemo(() => moleculeIds.slice(0, 500), [moleculeIds]);
  const eventSuggestions = useMemo(() => eventIds.slice(0, 500), [eventIds]);
  const [draft, setDraft] = useState<EvidenceSearchFilters>({});
  const [filters, setFilters] = useState<EvidenceSearchFilters>({});
  const [cursor, setCursor] = useState(0);
  const search = useQuery({
    queryKey: ['evidence-search', jobId, filters.moleculeId ?? '', filters.eventId ?? '', filters.state ?? '', cursor],
    queryFn: ({ signal }) => searchEvidence(jobId, filters, cursor, REMOTE_SEARCH_LIMIT, signal),
    staleTime: Number.POSITIVE_INFINITY,
    gcTime: 5 * 60 * 1000
  });
  const rows = search.data?.rows ?? [];
  const total = search.data?.total_matches ?? 0;
  const firstShown = rows.length === 0 ? 0 : cursor + 1;
  const lastShown = cursor + rows.length;
  const applyFilters = () => {
    setCursor(0);
    setFilters({
      moleculeId: normalizeExactEvidenceFilter(draft.moleculeId),
      eventId: normalizeExactEvidenceFilter(draft.eventId),
      state: draft.state || undefined
    });
  };
  return (
    <section className="overflow-hidden rounded-md border border-line" aria-labelledby="global-evidence-search-title">
      <div className="border-b border-line bg-panel2 px-3 py-3">
        <h3 id="global-evidence-search-title" className="text-sm font-semibold">Global evidence search</h3>
        <p className="mt-0.5 text-xs text-muted">
          Exact filters are evaluated over every immutable page by the server index; empty filters browse all stored observations. Suggestions are bounded to 500 IDs, but any exact ID is accepted.
        </p>
        <form
          className="mt-3 grid gap-2 md:grid-cols-[1fr_1fr_0.7fr_auto_auto]"
          onSubmit={(event) => { event.preventDefault(); applyFilters(); }}
        >
          <label className="grid gap-1 text-xs text-muted">
            Molecule ID
            <input
              type="text"
              list="evidence-molecule-ids"
              value={draft.moleculeId ?? ''}
              onChange={(event) => setDraft((current) => ({ ...current, moleculeId: event.target.value }))}
              className="rounded border border-line bg-shell px-2 py-1.5 text-sm text-text"
              placeholder="exact molecule ID"
            />
          </label>
          <datalist id="evidence-molecule-ids">
            {moleculeSuggestions.map((id) => <option key={id} value={id} />)}
          </datalist>
          <label className="grid gap-1 text-xs text-muted">
            Event ID
            <input
              type="text"
              list="evidence-event-ids"
              value={draft.eventId ?? ''}
              onChange={(event) => setDraft((current) => ({ ...current, eventId: event.target.value }))}
              className="rounded border border-line bg-shell px-2 py-1.5 font-mono text-sm text-text"
              placeholder="exact event ID"
            />
          </label>
          <datalist id="evidence-event-ids">
            {eventSuggestions.map((id) => <option key={id} value={id} />)}
          </datalist>
          <label className="grid gap-1 text-xs text-muted">
            Stored state
            <select
              value={draft.state ?? ''}
              onChange={(event) => setDraft((current) => ({ ...current, state: event.target.value || undefined }))}
              className="rounded border border-line bg-shell px-2 py-1.5 text-sm text-text"
            >
              <option value="">all states</option>
              {['ALTERNATE', 'REFERENCE', 'EVENT_ABSENT', 'LOW_QUALITY', 'CONFLICT'].map((state) => (
                <option key={state} value={state}>{state}</option>
              ))}
            </select>
          </label>
          <button type="submit" className="self-end rounded border border-aqua px-3 py-1.5 text-sm text-aqua">
            Apply
          </button>
          <button
            type="button"
            className="self-end rounded border border-line px-3 py-1.5 text-sm text-muted"
            onClick={() => { setDraft({}); setFilters({}); setCursor(0); }}
          >
            Clear
          </button>
        </form>
      </div>
      <div className="flex items-center justify-between gap-3 border-b border-line px-3 py-2 text-xs text-muted" aria-live="polite">
        <span>{search.isFetching ? 'Searching…' : `${firstShown}–${lastShown} of ${total} exact matches`}</span>
        <span className="flex gap-2">
          <button
            type="button"
            className="rounded border border-line px-2 py-1 disabled:opacity-40"
            disabled={cursor === 0 || search.isFetching}
            onClick={() => setCursor(Math.max(0, cursor - REMOTE_SEARCH_LIMIT))}
          >
            Previous matches
          </button>
          <button
            type="button"
            className="rounded border border-line px-2 py-1 disabled:opacity-40"
            disabled={search.data?.next_cursor == null || search.isFetching}
            onClick={() => setCursor(search.data?.next_cursor ?? cursor)}
          >
            Next matches
          </button>
        </span>
      </div>
      {search.isLoading ? <EmptyRow text="Building the first bounded result page…" /> : search.isError ? (
        <div role="alert" className="px-3 py-4 text-sm text-red-300">{(search.error as Error).message}</div>
      ) : rows.length === 0 ? <EmptyRow text="No stored observation matches the exact filters; this does not imply a reference state" /> : (
        <>
          <div className="grid grid-cols-[0.55fr_1fr_1fr_0.7fr_0.9fr_1.1fr] gap-2 border-b border-line px-3 py-2 text-xs font-semibold uppercase tracking-normal text-muted">
            <div>Page:row</div><div>Molecule</div><div>Event</div><div>State</div><div>Alignment</div><div>Evidence</div>
          </div>
          <FixedSizeList
            height={Math.min(320, Math.max(42, rows.length * 42))}
            width="100%"
            itemCount={rows.length}
            itemSize={42}
            itemData={rows}
            itemKey={(index, items) => items[index].id}
          >
            {EvidenceSearchResultRow}
          </FixedSizeList>
        </>
      )}
    </section>
  );
}

function EvidenceSearchResultRow({ index, style, data }: ListChildComponentProps<EvidenceSearchRow[]>) {
  const observation = data[index];
  return (
    <div style={style} className="grid grid-cols-[0.55fr_1fr_1fr_0.7fr_0.9fr_1.1fr] items-center gap-2 border-b border-line px-3 text-sm">
      <div className="tabular-nums text-muted">{observation.page_index}:{observation.row_index}</div>
      <div className="truncate font-medium" title={observation.molecule_id}>{observation.molecule_id}</div>
      <div className="truncate" title={observation.event_id}>{observation.event_id}</div>
      <div>{observation.state}</div>
      <div className="truncate" title={observation.alignment_id}>{observation.alignment_id}</div>
      <div className="truncate text-muted" title={observation.evidence_source}>{observation.evidence_source}</div>
    </div>
  );
}

function RemoteEvidenceBrowser({
  page,
  pageIndex,
  pageCount,
  loading,
  error,
  onPageChange
}: {
  page?: NonNullable<MitoAnalysisData['observation_pages']>[number];
  pageIndex: number;
  pageCount: number;
  loading: boolean;
  error: Error | null;
  onPageChange: (page: number) => void;
}) {
  const rows = useMemo(() => page ? materializeObservationPage(page) : [], [page]);
  const lastPage = Math.max(0, pageCount - 1);
  return (
    <section className="overflow-hidden rounded-md border border-line" aria-labelledby="remote-evidence-title">
      <div className="flex flex-wrap items-center justify-between gap-3 border-b border-line bg-panel2 px-3 py-2">
        <div>
          <h3 id="remote-evidence-title" className="text-sm font-semibold">Remote evidence pages</h3>
          <p className="mt-0.5 text-xs text-muted">Only the open immutable page is allocated and rendered in the browser.</p>
        </div>
        <div className="flex items-center gap-2 text-xs" aria-label="Evidence page navigation">
          <button type="button" className="rounded border border-line px-2 py-1 disabled:opacity-40" aria-label="First evidence page" disabled={pageIndex === 0 || pageCount === 0} onClick={() => onPageChange(0)}>First</button>
          <button type="button" className="rounded border border-line px-2 py-1 disabled:opacity-40" aria-label="Previous evidence page" disabled={pageIndex === 0 || pageCount === 0} onClick={() => onPageChange(pageIndex - 1)}>Prev</button>
          <label className="flex items-center gap-1 text-muted">
            Page
            <input
              className="w-16 rounded border border-line bg-shell px-2 py-1 text-right text-text"
              type="number"
              min={1}
              max={Math.max(1, pageCount)}
              value={pageCount === 0 ? 0 : pageIndex + 1}
              disabled={pageCount === 0}
              aria-label="Evidence page number"
              onChange={(event) => {
                const requested = Number.parseInt(event.target.value, 10);
                if (Number.isFinite(requested)) onPageChange(Math.min(lastPage, Math.max(0, requested - 1)));
              }}
            />
            / {pageCount}
          </label>
          <button type="button" className="rounded border border-line px-2 py-1 disabled:opacity-40" aria-label="Next evidence page" disabled={pageIndex >= lastPage || pageCount === 0} onClick={() => onPageChange(pageIndex + 1)}>Next</button>
          <button type="button" className="rounded border border-line px-2 py-1 disabled:opacity-40" aria-label="Last evidence page" disabled={pageIndex >= lastPage || pageCount === 0} onClick={() => onPageChange(lastPage)}>Last</button>
        </div>
      </div>
      <div className="grid grid-cols-[1fr_1fr_0.7fr_0.9fr_1.2fr] gap-2 border-b border-line px-3 py-2 text-xs font-semibold uppercase tracking-normal text-muted">
        <div>Molecule</div><div>Event</div><div>State</div><div>Alignment</div><div>Evidence</div>
      </div>
      {loading ? <EmptyRow text="Loading immutable evidence page…" /> : error ? (
        <div role="alert" className="px-3 py-4 text-sm text-red-300">{error.message}</div>
      ) : rows.length === 0 ? <EmptyRow text="Evidence page is empty" /> : (
        <FixedSizeList
          height={Math.min(320, Math.max(42, rows.length * 42))}
          width="100%"
          itemCount={rows.length}
          itemSize={42}
          itemData={rows}
          itemKey={(index, items) => items[index].id}
        >
          {RemoteEvidenceRow}
        </FixedSizeList>
      )}
    </section>
  );
}

function RemoteEvidenceRow({ index, style, data }: ListChildComponentProps<EvidenceObservation[]>) {
  const observation = data[index];
  return (
    <div style={style} className="grid grid-cols-[1fr_1fr_0.7fr_0.9fr_1.2fr] items-center gap-2 border-b border-line px-3 text-sm">
      <div className="truncate font-medium" title={observation.molecule_id}>{observation.molecule_id}</div>
      <div className="truncate" title={observation.event_id}>{observation.event_id}</div>
      <div>{observation.state}</div>
      <div className="truncate" title={observation.alignment_id}>{observation.alignment_id}</div>
      <div className="truncate text-muted" title={observation.evidence_source}>{observation.evidence_source}</div>
    </div>
  );
}

function materializeObservationPage(
  page: NonNullable<MitoAnalysisData['observation_pages']>[number]
): EvidenceObservation[] {
  return Array.from({ length: page.count }, (_, index) => ({
    id: `observation:${page.offset + index}`,
    molecule_id: page.columns.molecule_id[index],
    event_id: page.columns.event_id[index],
    alignment_id: page.columns.alignment_id[index],
    state: page.columns.state[index],
    observed_allele: page.columns.observed_allele[index],
    base_quality: page.columns.base_quality[index],
    mapping_quality: page.columns.mapping_quality[index],
    strand: page.columns.strand[index],
    evidence_source: page.columns.evidence_source[index],
    read_position: page.columns.read_position[index]
  }));
}

function RearrangementsPanel({ data }: { data: MitoAnalysisData }) {
  const paths = data.complex_events ?? [];
  const eventsById = new Map((data.events ?? []).map((event) => [event.id, event]));
  const [selectedPathId, setSelectedPathId] = useState<string | undefined>();
  const selectedPath = paths.find((path) => path.id === selectedPathId);
  const selectedEvent = selectedPath
    ? eventsById.get(selectedPath.event_id ?? selectedPath.id)
    : undefined;
  const setSelectedEventId = useMitoStore((state) => state.setSelectedEventId);
  const setSelectedMoleculeId = useMitoStore((state) => state.setSelectedMoleculeId);
  const setSelectedPhaseId = useMitoStore((state) => state.setSelectedPhaseId);
  const setSelectedSvId = useMitoStore((state) => state.setSelectedSvId);
  const supportingMoleculeIds = selectedPath
    ? selectedEvent?.supporting_molecule_ids ?? selectedPath.supporting_reads
    : [];
  const supportingMoleculeSet = new Set(supportingMoleculeIds);
  const supportingAlignments = (data.alignments ?? []).filter((alignment) =>
    supportingMoleculeSet.has(alignment.molecule_id)
  );
  const inspectEvent = (eventId: string) => {
    setSelectedMoleculeId(undefined);
    setSelectedPhaseId(undefined);
    setSelectedEventId(eventId);
  };
  const inspectMolecule = (moleculeId: string) => {
    setSelectedEventId(undefined);
    setSelectedPhaseId(undefined);
    setSelectedMoleculeId(moleculeId);
  };

  return (
    <div className="grid gap-4">
      <div className="grid gap-3 sm:grid-cols-3">
        <Metric label="Structural events" value={data.svs.length.toString()} />
        <Metric label="Multi-junction paths" value={paths.length.toString()} />
        <Metric
          label="Path-support molecules"
          value={new Set(paths.flatMap((path) => eventsById.get(path.event_id ?? path.id)?.supporting_molecule_ids ?? path.supporting_reads)).size.toString()}
        />
      </div>

      <section className="overflow-hidden rounded-md border border-line" aria-labelledby="complex-paths-title">
        <div className="border-b border-line bg-panel2 px-3 py-2">
          <h3 id="complex-paths-title" className="text-sm font-semibold">Observed multi-junction paths</h3>
          <p className="mt-0.5 text-xs text-muted">
            Ordered split-alignment evidence only. Missing junctions and circular closure are never inferred.
          </p>
        </div>
        <div className="grid grid-cols-[1.6fr_0.55fr_0.55fr_0.7fr_0.75fr] gap-2 border-b border-line px-3 py-2 text-xs font-semibold uppercase tracking-normal text-muted">
          <div>Canonical path</div><div>Junctions</div><div>Segments</div><div>Molecules</div><div>Inspect</div>
        </div>
        {paths.length === 0 ? <EmptyRow text="No molecule contains two or more observed split-alignment junctions" /> : paths.map((path) => {
          const event = eventsById.get(path.event_id ?? path.id);
          const support = event?.supporting_molecule_ids ?? path.supporting_reads;
          const selected = selectedPath?.id === path.id;
          return (
            <button
              key={path.id}
              type="button"
              aria-pressed={selected}
              onClick={() => {
                setSelectedPathId(path.id);
                inspectEvent(path.event_id ?? path.id);
              }}
              className="grid w-full grid-cols-[1.6fr_0.55fr_0.55fr_0.7fr_0.75fr] gap-2 border-b border-line px-3 py-2 text-left text-sm last:border-0 hover:bg-panel2 aria-pressed:bg-panel2"
            >
              <span className="truncate font-mono text-xs" title={path.id}>{path.id}</span>
              <span>{path.junction_count}</span>
              <span>{path.segment_count}</span>
              <span>{support.length}</span>
              <span className="text-aqua">{selected ? 'selected' : 'trace path'}</span>
            </button>
          );
        })}
      </section>

      {selectedPath && (
        <section className="overflow-hidden rounded-md border border-aqua/50" aria-labelledby="selected-path-title">
          <div className="flex flex-wrap items-start justify-between gap-3 border-b border-line bg-panel2 px-3 py-3">
            <div className="min-w-0">
              <h3 id="selected-path-title" className="font-semibold">Selected rearrangement path</h3>
              <p className="mt-1 break-all font-mono text-xs text-aqua">{selectedPath.id}</p>
            </div>
            <button
              type="button"
              className="text-sm text-aqua hover:underline"
              onClick={() => { setSelectedPathId(undefined); setSelectedEventId(undefined); }}
            >
              Clear selection
            </button>
          </div>
          <div className="grid gap-4 p-3 lg:grid-cols-[minmax(0,1fr)_minmax(320px,0.8fr)]">
            <div>
              <div className="text-xs font-semibold uppercase tracking-normal text-muted">Ordered junction evidence</div>
              <ol className="mt-2 grid gap-2" aria-label="Junction traversal order">
                {selectedPath.junction_ids.map((junctionId, index) => {
                  const eventId = `sv:${junctionId}`;
                  const sv = data.svs.find((candidate) => candidate.id === junctionId);
                  return (
                    <li key={`${junctionId}:${index}`} className="grid grid-cols-[auto_1fr_auto] items-center gap-3 rounded border border-line bg-shell px-3 py-2">
                      <span className="grid h-7 w-7 place-items-center rounded-full border border-aqua/50 text-xs text-aqua">{index + 1}</span>
                      <span className="min-w-0">
                        <button
                          type="button"
                          className="block max-w-full truncate text-left font-mono text-xs hover:text-aqua"
                          title={junctionId}
                          onClick={() => {
                            setSelectedSvId(junctionId);
                            inspectEvent(sv?.event_id ?? eventId);
                          }}
                        >
                          {junctionId}
                        </button>
                        <span className="mt-0.5 block text-xs text-muted">
                          {sv ? `${sv.type} · ${sv.start}–${sv.end} · ${sv.length} bp` : 'component projection unavailable'}
                        </span>
                      </span>
                      <span className="font-mono text-xs text-amber" title="Observed strand transition">
                        {selectedPath.junction_orientations[index] ?? 'unknown'}
                      </span>
                    </li>
                  );
                })}
              </ol>
              <div className="mt-3 rounded border border-amber/50 bg-amber/10 px-3 py-2 text-xs text-amber">
                Canonicalization: {selectedPath.canonicalization}. This is an observed open traversal; rotation equivalence and molecular closure are not claimed.
              </div>
            </div>
            <div className="grid content-start gap-3">
              <TraceList
                label="Supporting molecules"
                values={supportingMoleculeIds}
                empty="No molecule reference was retained"
                onSelect={inspectMolecule}
              />
              <div className="rounded border border-line bg-shell p-3">
                <div className="text-xs font-semibold uppercase tracking-normal text-muted">Evidence contract</div>
                <dl className="mt-2 grid grid-cols-[auto_1fr] gap-x-3 gap-y-1 text-xs">
                  <dt className="text-muted">Event</dt><dd className="break-all font-mono">{selectedPath.event_id ?? selectedPath.id}</dd>
                  <dt className="text-muted">Assessability</dt><dd>{selectedEvent?.assessability ?? 'legacy projection'}</dd>
                  <dt className="text-muted">Negative evidence</dt><dd>{selectedEvent?.negative_evidence_rule ?? 'not encoded'}</dd>
                  <dt className="text-muted">Topology</dt><dd>observed partial path; closure not inferred</dd>
                </dl>
              </div>
            </div>
          </div>
          <AlignmentTraceTable alignments={supportingAlignments} />
        </section>
      )}

      <section className="overflow-hidden rounded-md border border-line" aria-labelledby="simple-sv-title">
        <div className="border-b border-line bg-panel2 px-3 py-2">
          <h3 id="simple-sv-title" className="text-sm font-semibold">Canonical structural events</h3>
        </div>
        <div className="grid grid-cols-[1.2fr_0.7fr_0.7fr_0.6fr_0.65fr_1fr] gap-2 border-b border-line px-3 py-2 text-xs font-semibold uppercase tracking-normal text-muted">
          <div>Event</div><div>Type</div><div>Coordinates</div><div>Length</div><div>Molecules</div><div>Evidence / orientation</div>
        </div>
        {data.svs.length === 0 ? <EmptyRow text="No canonical structural event" /> : data.svs.map((sv) => {
          const event = eventsById.get(sv.event_id ?? `sv:${sv.id}`);
          return (
            <button
              key={sv.id}
              type="button"
              onClick={() => { setSelectedSvId(sv.id); inspectEvent(sv.event_id ?? `sv:${sv.id}`); }}
              className="grid w-full grid-cols-[1.2fr_0.7fr_0.7fr_0.6fr_0.65fr_1fr] gap-2 border-b border-line px-3 py-2 text-left text-sm last:border-0 hover:bg-panel2"
            >
              <span className="truncate font-mono text-xs" title={sv.id}>{sv.id}</span>
              <span>{sv.type}</span>
              <span>{sv.start}–{sv.end}</span>
              <span>{sv.length}</span>
              <span>{event?.supporting_molecule_ids.length ?? sv.supporting_reads.length}</span>
              <span className="truncate text-muted" title={[...(sv.evidence_sources ?? []), ...(sv.orientations ?? [])].join(', ')}>
                {(sv.evidence_sources ?? [sv.evidence_source ?? 'unknown']).join('+')} · {(sv.orientations ?? [sv.orientation ?? 'n/a']).join(', ')}
              </span>
            </button>
          );
        })}
      </section>
    </div>
  );
}

function TraceList({
  label,
  values,
  empty,
  onSelect
}: {
  label: string;
  values: string[];
  empty: string;
  onSelect: (value: string) => void;
}) {
  return (
    <div className="rounded border border-line bg-shell p-3">
      <div className="text-xs font-semibold uppercase tracking-normal text-muted">{label}</div>
      <div className="mt-2 max-h-36 overflow-y-auto">
        {values.length === 0 ? <span className="text-xs text-muted">{empty}</span> : values.map((value) => (
          <button
            key={value}
            type="button"
            className="block w-full truncate rounded px-1 py-1 text-left font-mono text-xs hover:bg-panel2 hover:text-aqua"
            title={value}
            onClick={() => onSelect(value)}
          >
            {value}
          </button>
        ))}
      </div>
    </div>
  );
}

function AlignmentTraceTable({ alignments }: { alignments: AlignmentFragment[] }) {
  return (
    <div className="border-t border-line">
      <div className="px-3 py-2 text-xs font-semibold uppercase tracking-normal text-muted">Source alignment fragments</div>
      <div className="grid grid-cols-[0.8fr_0.9fr_0.55fr_0.6fr_0.65fr_0.8fr_1.4fr] gap-2 border-t border-line bg-panel2 px-3 py-2 text-xs font-semibold uppercase tracking-normal text-muted">
        <div>Alignment</div><div>Molecule</div><div>Role</div><div>Strand</div><div>MAPQ</div><div>Start / CIGAR</div><div>SA provenance</div>
      </div>
      {alignments.length === 0 ? <EmptyRow text="Alignment-fragment provenance is unavailable in this result schema" /> : alignments.map((alignment) => (
        <div key={alignment.id} className="grid grid-cols-[0.8fr_0.9fr_0.55fr_0.6fr_0.65fr_0.8fr_1.4fr] gap-2 border-t border-line px-3 py-2 text-xs">
          <div className="truncate font-mono" title={alignment.id}>{alignment.id}</div>
          <div className="truncate" title={alignment.molecule_id}>{alignment.molecule_id}</div>
          <div>{alignment.role}</div>
          <div>{alignment.strand}</div>
          <div>{alignment.mapping_quality}</div>
          <div className="truncate" title={`${alignment.reference_name}:${alignment.reference_start} ${alignment.cigar}`}>
            {alignment.reference_start} · {alignment.cigar}
          </div>
          <div className="truncate font-mono text-muted" title={alignment.aux_tags?.SA ?? 'no SA tag'}>
            {alignment.aux_tags?.SA ?? 'no SA tag retained'}
          </div>
        </div>
      ))}
    </div>
  );
}

function HaplogroupPanel({ data }: { data: MitoAnalysisData }) {
  return (
    <div className="overflow-hidden rounded-md border border-line">
      <div className="grid grid-cols-[0.7fr_0.8fr_0.55fr_0.8fr_1.4fr] gap-2 border-b border-line bg-panel2 px-3 py-2 text-xs font-semibold uppercase tracking-normal text-muted">
        <div>Cluster</div>
        <div>Haplogroup</div>
        <div>Quality</div>
        <div>Molecules</div>
        <div>Alternatives / molecule signature</div>
      </div>
      {data.clusters.length === 0 ? (
        <EmptyRow text="No clusters" />
      ) : (
        data.clusters.map((cluster) => (
          <div
            key={cluster.id}
            className="grid grid-cols-[0.7fr_0.8fr_0.55fr_0.8fr_1.4fr] gap-2 border-b border-line px-3 py-2 text-sm last:border-0"
          >
            <div className="font-semibold">{cluster.label}</div>
            <div className={cluster.haplogroup === 'unassigned' ? 'text-muted' : 'text-aqua'}>
              {cluster.haplogroup ?? 'unassigned'}
              {cluster.haplogroup_assignment?.contamination_warning && (
                <span className="ml-2 text-amber" title="Competing macrohaplogroup assignments">
                  mixed?
                </span>
              )}
            </div>
            <div>{cluster.haplogroup_assignment ? `${cluster.haplogroup_assignment.quality.toFixed(1)}%` : '—'}</div>
            <div>{cluster.size}</div>
            <div className="truncate text-muted">
              {cluster.haplogroup_assignment?.candidates.slice(1).map((candidate) =>
                `${candidate.name} ${candidate.score.toFixed(1)}%`
              ).join(', ') || 'no ranked alternatives'}
              {' | '}
              {cluster.sv_signature.length === 0 ? 'no SVs' : cluster.sv_signature.map((sv) => `${sv.sv_id} n=${sv.support}`).join(', ')}
              {(cluster.complex_event_signature?.length ?? 0) > 0
                ? ` | complex paths: ${cluster.complex_event_signature?.map((event) => `n=${event.support}`).join(', ')}`
                : ''}
              {cluster.haplogroup_assignment
                ? ` | markers: ${cluster.haplogroup_assignment.observed_markers?.length ?? 0}, callable bp: ${cluster.haplogroup_assignment.callable_ranges?.reduce((total, range) => total + range.end - range.start + 1, 0) ?? 'unknown'}`
                : ''}
            </div>
          </div>
        ))
      )}
    </div>
  );
}

function ClinicalPanel({ variants }: { variants: VariantSummary[] }) {
  return (
    <div className="overflow-hidden rounded-md border border-line">
      <div className="grid grid-cols-[1fr_0.55fr_0.65fr_1fr_0.75fr_0.6fr] gap-2 border-b border-line bg-panel2 px-3 py-2 text-xs font-semibold uppercase tracking-normal text-muted">
        <div>Variant</div>
        <div>Support</div>
        <div>Frequency</div>
        <div>Clinical summary</div>
        <div>Conflict</div>
        <div>Records</div>
      </div>
      {variants.length === 0 ? (
        <EmptyRow text="No variants" />
      ) : (
        variants.map((variant) => {
          const annotation = variant.annotation;
          const assertions = annotation?.assertions ?? [];
          const conflicting = annotation?.conflict_status === 'conflicting';
          return (
            <div key={variant.key} className="border-b border-line last:border-0">
              <div className="grid grid-cols-[1fr_0.55fr_0.65fr_1fr_0.75fr_0.6fr] gap-2 px-3 py-2 text-sm">
                <div className="font-semibold">{variant.label}</div>
                <div>{variant.support}</div>
                <div>{formatPercent(variant.frequency)}</div>
                <div className={annotation ? 'text-amber' : 'text-muted'}>
                  {annotation?.consensus_significance ?? annotation?.pathogenicity ?? 'not annotated'}
                </div>
                <div className={conflicting ? 'font-semibold text-red-300' : annotation ? 'text-text' : 'text-muted'}>
                  {annotation?.conflict_status ?? 'not assessed'}
                </div>
                <div>{assertions.length || '—'}</div>
              </div>
              {annotation && (
                <details className="border-t border-line/70 bg-shell/40 px-3 py-2 text-xs text-muted">
                  <summary className="cursor-pointer font-semibold text-text">
                    Source assertions and provenance
                  </summary>
                  <div className="mt-2 grid gap-2">
                    {assertions.length === 0 ? (
                      <div>Legacy summary only; assertion-level provenance is unavailable.</div>
                    ) : (
                      assertions.map((assertion, index) => (
                        <div
                          key={`${assertion.source}:${assertion.assertion_id}:${index}`}
                          className="rounded border border-line bg-panel2 p-2"
                        >
                          <div className="flex flex-wrap items-center gap-x-2 gap-y-1 text-text">
                            <span className="font-semibold">{assertion.source}</span>
                            <span>{assertion.clinical_significance || 'significance not provided'}</span>
                            <span className="text-muted">({assertion.normalized_significance})</span>
                          </div>
                          <div className="mt-1">
                            {assertion.disease || 'disease not provided'}
                            {assertion.review_status ? ` | review: ${assertion.review_status}` : ''}
                            {assertion.assertion_date ? ` | asserted: ${assertion.assertion_date}` : ''}
                            {assertion.allele_id ? ` | allele ${assertion.allele_id}` : ''}
                          </div>
                          <div className="mt-1 break-words">
                            snapshot {assertion.resource_version || 'not provided'} / {assertion.retrieved_at || 'retrieval unknown'}
                            {assertion.references.length ? ` | ${assertion.references.join(', ')}` : ''}
                            {assertion.source_url && (
                              <>
                                {' | '}
                                <a
                                  className="text-aqua underline underline-offset-2"
                                  href={assertion.source_url}
                                  target="_blank"
                                  rel="noreferrer"
                                >
                                  source record
                                </a>
                              </>
                            )}
                          </div>
                        </div>
                      ))
                    )}
                  </div>
                </details>
              )}
            </div>
          );
        })
      )}
    </div>
  );
}

function CoveragePanel({
  data,
  metrics
}: {
  data: MitoAnalysisData;
  metrics: ReturnType<typeof summarizeCoverage>;
}) {
  const preview = data.coverage.slice(0, 80);
  return (
    <div className="grid gap-4 lg:grid-cols-[320px_minmax(0,1fr)]">
      <div className="grid gap-3 sm:grid-cols-3 lg:grid-cols-1">
        <Metric label="Mean depth" value={metrics.meanDepth.toFixed(2)} />
        <Metric label="Bins >20x" value={formatPercent(metrics.pctGt20 / 100)} />
        <Metric label="Max depth" value={metrics.maxDepth.toString()} />
      </div>
      <div className="rounded-md border border-line bg-panel2 p-3">
        <div className="mb-2 text-xs font-semibold uppercase tracking-normal text-muted">Coverage histogram</div>
        <div className="flex h-28 items-end gap-px overflow-hidden">
          {preview.map((bin) => (
            <div
              key={`${bin.start}-${bin.end}`}
              className="min-w-[3px] flex-1 rounded-t-sm bg-sky"
              title={`${bin.start}-${bin.end}: ${bin.depth}`}
              style={{ height: `${Math.max(4, (bin.depth / Math.max(1, metrics.maxDepth)) * 104)}px` }}
            />
          ))}
        </div>
      </div>
    </div>
  );
}

function ProteinPanel({ variants }: { variants: VariantSummary[] }) {
  const proteinVariants = variants.filter((variant) => variant.protein || variant.structure?.structure_id);
  const [selectedKey, setSelectedKey] = useState<string>();
  const selectedEventId = useMitoStore((state) => state.selectedEventId);
  const setSelectedEventId = useMitoStore((state) => state.setSelectedEventId);
  const selected = proteinVariants.find((variant) => variant.key === selectedKey) ?? proteinVariants[0];

  useEffect(() => {
    if (selectedEventId && proteinVariants.some((variant) => variant.key === selectedEventId)) {
      setSelectedKey(selectedEventId);
    } else if (!selectedKey && proteinVariants[0]) {
      setSelectedKey(proteinVariants[0].key);
    } else if (selectedKey && !proteinVariants.some((variant) => variant.key === selectedKey)) {
      setSelectedKey(proteinVariants[0]?.key);
    }
  }, [proteinVariants, selectedEventId, selectedKey]);

  return (
    <div className="grid gap-4 lg:grid-cols-[320px_minmax(0,1fr)]">
      <div className="grid content-start gap-3">
        <Metric label="Mapped residues" value={proteinVariants.length.toString()} />
        <Metric label="Viewer" value={selected?.structure?.structure_id ?? 'NGL-ready'} />
        <Metric label="Example" value={selected ? selected.residue ?? selected.label : 'no protein variant'} />
        <div className="overflow-hidden rounded-md border border-line bg-panel2">
          <div className="border-b border-line px-3 py-2 text-xs font-semibold uppercase tracking-normal text-muted">
            Protein variants
          </div>
          <div className="max-h-52 overflow-auto scrollbar-thin">
            {proteinVariants.length === 0 ? (
              <EmptyRow text="No protein-mapped variants" />
            ) : (
              proteinVariants.map((variant) => (
                <button
                  key={variant.key}
                  type="button"
                  onClick={() => {
                    setSelectedKey(variant.key);
                    setSelectedEventId(variant.key);
                  }}
                  className={[
                    'w-full border-b px-3 py-2 text-left text-sm last:border-0',
                    selected?.key === variant.key ? 'border-aqua bg-teal-300/10' : 'border-line hover:bg-panel'
                  ].join(' ')}
                >
                  <div className="font-semibold">{variant.protein ?? variant.gene}</div>
                  <div className="mt-1 truncate text-muted">{variant.residue ?? variant.label}</div>
                </button>
              ))
            )}
          </div>
        </div>
      </div>
      <div className="rounded-md border border-line bg-panel2 p-4">
        {selected ? (
          <div className="grid gap-4">
            <NglProteinViewer variant={toProteinViewerVariant(selected)} />
            <dl className="grid gap-3 sm:grid-cols-2 xl:grid-cols-4">
              <Metric label="Variant" value={selected.label} />
              <Metric label="Residue" value={selected.residue ?? 'mapped'} />
              <Metric label="Structure" value={selected.structure?.structure_id ?? 'local model'} />
              <Metric label="Frequency" value={formatPercent(selected.frequency)} />
            </dl>
            <MoleculeDistribution variant={selected} />
          </div>
        ) : (
          <div className="text-sm text-muted">
            Non-synonymous consequence calls and local structure mappings are required before residue highlights are populated.
          </div>
        )}
      </div>
    </div>
  );
}

function MoleculeDistribution({ variant }: { variant: VariantSummary }) {
  const clusters = new Map<number, number>();
  for (const molecule of variant.molecules) {
    clusters.set(molecule.clusterId, (clusters.get(molecule.clusterId) ?? 0) + 1);
  }
  const groups = [...clusters.entries()].sort(([left], [right]) => left - right);
  const maximum = Math.max(1, ...groups.map(([, count]) => count));

  return (
    <div className="rounded-md border border-line bg-panel p-3">
      <div className="flex flex-wrap items-center justify-between gap-2">
        <div>
          <div className="text-xs font-semibold uppercase tracking-normal text-muted">Molecule-level support</div>
          <div className="mt-1 text-sm text-text">
            {variant.molecules.length} observed molecules across {groups.length} read cluster{groups.length === 1 ? '' : 's'}
          </div>
        </div>
        <div className="text-xs text-muted">
          {variant.callableDepth == null
            ? 'Observed fraction; locus-callable depth unavailable'
            : `HF ${formatPercent(variant.frequency)} | DP ${variant.callableDepth} | 95% CI ${formatPercent(
                variant.ci95Low ?? 0
              )}-${formatPercent(variant.ci95High ?? 0)}`}
        </div>
      </div>
      {groups.length > 0 && (
        <div className="mt-3 grid gap-2 sm:grid-cols-2 xl:grid-cols-3">
          {groups.map(([clusterId, count]) => (
            <div key={clusterId} className="rounded border border-line bg-panel2 px-2.5 py-2 text-xs">
              <div className="flex items-center justify-between gap-2">
                <span className="font-semibold text-text">{clusterId < 0 ? 'Outliers' : `Cluster ${clusterId + 1}`}</span>
                <span className="text-aqua">n={count}</span>
              </div>
              <div className="mt-2 h-1.5 overflow-hidden rounded-full bg-shell">
                <div className="h-full rounded-full bg-aqua" style={{ width: `${(count / maximum) * 100}%` }} />
              </div>
            </div>
          ))}
        </div>
      )}
      <details className="mt-3 text-xs text-muted">
        <summary className="cursor-pointer select-none font-semibold text-text">Supporting molecule IDs</summary>
        <div className="mt-2 max-h-28 overflow-auto rounded border border-line bg-shell p-2 font-mono">
          {variant.molecules.map((molecule) => (
            <div key={`${molecule.clusterId}:${molecule.id}`} className="truncate" title={molecule.id}>
              {molecule.id}
            </div>
          ))}
        </div>
      </details>
      {(variant.strandSupport || variant.alleleQuality || variant.readPosition) && (
        <div className="mt-3 grid gap-2 sm:grid-cols-2 xl:grid-cols-4">
          <Metric
            label="Alt strand F/R"
            value={variant.strandSupport ? `${variant.strandSupport.alt_forward}/${variant.strandSupport.alt_reverse}` : '—'}
          />
          <Metric
            label="Strand delta"
            value={variant.strandBiasDelta == null ? 'not estimable' : variant.strandBiasDelta.toFixed(3)}
          />
          <Metric
            label="Alt mean Phred"
            value={variant.alleleQuality?.alternate.mean_phred?.toFixed(1) ?? 'not estimable'}
          />
          <Metric
            label="Read-position delta"
            value={variant.readPosition?.bias_delta == null ? 'not estimable' : variant.readPosition.bias_delta.toFixed(3)}
          />
        </div>
      )}
    </div>
  );
}

function toProteinViewerVariant(variant: VariantSummary): ProteinViewerVariant {
  return {
    key: variant.key,
    label: variant.label,
    gene: variant.gene,
    protein: variant.protein,
    residue: variant.residue,
    frequency: variant.frequency,
    structure: variant.structure
  };
}

function Metric({ label, value }: { label: string; value: string }) {
  return (
    <div className="min-w-0 rounded-md border border-line bg-panel2 p-3">
      <div className="text-xs font-semibold uppercase tracking-normal text-muted">{label}</div>
      <div className="mt-1 break-words text-lg font-semibold leading-snug tracking-normal [overflow-wrap:anywhere]">
        {value}
      </div>
    </div>
  );
}

function EmptyRow({ text }: { text: string }) {
  return <div className="px-3 py-3 text-sm text-muted">{text}</div>;
}

function summarizeVariants(reads: ReadFeature[], data: MitoAnalysisData): VariantSummary[] {
  const denominator = Math.max(1, reads.length);
  const snps = new Map<string, VariantSummary>();
  const readClusterById = new Map(reads.map((read) => [read.id, read.cluster_id]));

  if (data.variants) {
    for (const snp of data.variants) {
      const key = snp.event_id ?? `snp:${snp.position}:${snp.ref}:${snp.alt}`;
      const type = snp.type ?? 'SNV';
      snps.set(key, {
        key,
        label: type === 'SNV'
          ? `${snp.position} ${snp.ref}>${snp.alt}`
          : `${type.replace('SMALL_', '').toLowerCase()} ${snp.position} ${snp.ref}>${snp.alt}`,
        type,
        position: snp.position,
        support: snp.alt_depth,
        frequency: snp.heteroplasmy,
        gene: snp.gene,
        consequence: snp.consequence,
        protein: snp.protein,
        residue: snp.residue,
        annotation: snp.annotation,
        structure: snp.structure,
        molecules: (snp.supporting_molecule_ids ?? snp.supporting_reads).map((id) => ({
          id,
          clusterId: readClusterById.get(id) ?? -1
        })),
        callableDepth: snp.callable_depth,
        ci95Low: snp.ci95_low,
        ci95High: snp.ci95_high,
        strandSupport: snp.strand_support,
        strandBiasDelta: snp.strand_bias_delta,
        alleleQuality: snp.allele_quality,
        readPosition: snp.read_position
      });
    }
  } else {
    for (const read of reads) {
      for (const snp of read.snps) {
      const key = `snp:${snp.position}:${snp.ref}:${snp.alt}`;
      const existing = snps.get(key);
      if (existing) {
        existing.support += 1;
        existing.frequency = existing.support / denominator;
        existing.molecules.push({ id: read.id, clusterId: read.cluster_id });
      } else {
        snps.set(key, {
          key,
          label: `${snp.position} ${snp.ref}>${snp.alt}`,
          type: 'SNP',
          position: snp.position,
          support: 1,
          frequency: 1 / denominator,
          gene: snp.gene,
          consequence: snp.consequence,
          protein: snp.protein,
          residue: snp.residue,
          annotation: snp.annotation,
          structure: snp.structure,
          molecules: [{ id: read.id, clusterId: read.cluster_id }]
        });
      }
    }
    }
  }

  const svs = data.svs.map((sv) => ({
    key: `sv:${sv.id}`,
    label: `${sv.type} ${sv.start}-${sv.end}`,
    type: sv.type,
    position: sv.start,
    support: sv.supporting_reads.length,
    frequency: sv.supporting_reads.length / denominator,
    annotation: sv.annotation,
    molecules: sv.supporting_reads.map((id) => ({
      id,
      clusterId: readClusterById.get(id) ?? -1
    }))
  }));

  return [...snps.values(), ...svs].sort((a, b) => a.position - b.position || a.label.localeCompare(b.label));
}

function summarizeCoverage(data: MitoAnalysisData) {
  if (data.coverage_metrics) {
    return {
      meanDepth: data.coverage_metrics.mean_depth,
      pctGt20: data.coverage_metrics.pct_sites_gt20x ?? data.coverage_metrics.pct_bins_gt20x ?? 0,
      maxDepth: data.coverage_metrics.max_depth
    };
  }

  const bins = data.coverage;
  const maxDepth = Math.max(0, ...bins.map((bin) => bin.depth));
  const meanDepth = bins.length === 0 ? 0 : bins.reduce((sum, bin) => sum + bin.depth, 0) / bins.length;
  const pctGt20 = bins.length === 0 ? 0 : (bins.filter((bin) => bin.depth > 20).length * 100) / bins.length;
  return { meanDepth, pctGt20, maxDepth };
}

function formatPercent(value: number): string {
  return `${(value * 100).toFixed(1)}%`;
}
