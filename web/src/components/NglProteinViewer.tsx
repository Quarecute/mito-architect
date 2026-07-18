import { Download, ExternalLink, Focus, Pause, Play, RotateCcw } from 'lucide-react';
import { useEffect, useMemo, useRef, useState } from 'react';
import type { ProteinStructureMapping } from '@mito-architect/visualization-lib';
import type { Stage, StructureComponent } from 'ngl';
import {
  normalizePdbId,
  rcsbCandidates,
  type StructureSourceCandidate
} from '../lib/proteinStructureSources';

export interface ProteinViewerVariant {
  key: string;
  label: string;
  gene?: string;
  protein?: string;
  residue?: string;
  frequency: number;
  structure?: ProteinStructureMapping;
}

type RepresentationMode = 'cartoon' | 'tube' | 'backbone' | 'surface';
type ColorMode = 'confidence' | 'secondary' | 'uniform';

interface NglProteinViewerProps {
  variant: ProteinViewerVariant;
}

interface LoadedModel {
  label: string;
  entryUrl: string;
  chain: string;
  coordinateDescription: string;
  confidenceSemantics: 'bfactor';
}

const MODES: Array<{ value: RepresentationMode; label: string }> = [
  { value: 'cartoon', label: 'Cartoon' },
  { value: 'tube', label: 'Tube' },
  { value: 'backbone', label: 'Backbone' },
  { value: 'surface', label: 'Surface' }
];

const COLORS: Array<{ value: ColorMode; label: string }> = [
  { value: 'confidence', label: 'Confidence' },
  { value: 'secondary', label: 'Secondary structure' },
  { value: 'uniform', label: 'Uniform' }
];

export default function NglProteinViewer({ variant }: NglProteinViewerProps) {
  const hostRef = useRef<HTMLDivElement>(null);
  const stageRef = useRef<Stage>();
  const componentRef = useRef<StructureComponent>();
  const loadGenerationRef = useRef(0);
  const [mode, setMode] = useState<RepresentationMode>('cartoon');
  const [colorMode, setColorMode] = useState<ColorMode>('confidence');
  const [showNeighborhood, setShowNeighborhood] = useState(true);
  const [showAssembly, setShowAssembly] = useState(false);
  const [spinning, setSpinning] = useState(false);
  const [status, setStatus] = useState('Initializing WebGL');
  const [stageReady, setStageReady] = useState(false);
  const [loadedModel, setLoadedModel] = useState<LoadedModel>();
  const residueIndex = variant.structure?.residue_index ?? parseResidueIndex(variant.residue);
  const mappedChain = sanitizeChain(variant.structure?.chain ?? 'A');
  const chain = loadedModel?.chain ?? mappedChain;
  const residueSelection = residueIndex ? `${residueIndex}:${chain}` : undefined;
  const contextSelection = residueSelection ? `within 5 of (${residueSelection})` : undefined;
  const structureId = variant.structure?.structure_id;

  const modelKey = useMemo(
    () => `${structureId ?? ''}|${mappedChain}|${residueIndex ?? ''}`,
    [mappedChain, residueIndex, structureId]
  );

  useEffect(() => {
    if (!hostRef.current) {
      return;
    }

    let disposed = false;
    let resizeObserver: ResizeObserver | undefined;
    let resize: (() => void) | undefined;

    import('ngl')
      .then(({ Stage }) => {
        if (disposed || !hostRef.current) {
          return;
        }

        const stage = new Stage(hostRef.current, {
          backgroundColor: '#08111f',
          cameraType: 'perspective',
          quality: 'high',
          sampleLevel: 1,
          tooltip: true
        });
        stageRef.current = stage;
        resizeObserver = new ResizeObserver(() => stage.handleResize());
        resizeObserver.observe(hostRef.current);
        resize = () => stage.handleResize();
        window.addEventListener('resize', resize);
        setStageReady(true);
      })
      .catch((error: unknown) => setStatus(errorMessage(error, 'NGL failed to load')));

    return () => {
      disposed = true;
      resizeObserver?.disconnect();
      if (resize) {
        window.removeEventListener('resize', resize);
      }
      stageRef.current?.dispose();
      stageRef.current = undefined;
      componentRef.current = undefined;
    };
  }, []);

  // Loading is intentionally isolated from representation and animation state. A model is
  // fetched once per variant instead of being downloaded again on every viewer interaction.
  useEffect(() => {
    const stage = stageRef.current;
    if (!stage || !stageReady) {
      return;
    }

    const generation = ++loadGenerationRef.current;
    setStatus('Resolving structure provenance');
    setLoadedModel(undefined);
    componentRef.current = undefined;
    stage.setSpin(false);
    stage.removeAllComponents();

    if (!structureId) {
      setStatus('No curated structure mapping is available for this variant');
      return;
    }

    loadProteinStructure(stage, structureId, mappedChain)
      .then(({ component, ...model }) => {
        if (generation !== loadGenerationRef.current) {
          stage.removeComponent(component);
          return;
        }

        componentRef.current = component;
        const loadedResidueSelection = residueIndex ? `${residueIndex}:${model.chain}` : undefined;
        const loadedContextSelection = loadedResidueSelection ? `within 5 of (${loadedResidueSelection})` : undefined;
        setLoadedModel(model);
        setStatus(
          loadedResidueSelection
            ? `Loaded; residue ${loadedResidueSelection} mapped`
            : 'Loaded; residue mapping unavailable'
        );
        component.autoView('all', 500);
        if (loadedResidueSelection) {
          window.setTimeout(() => {
            if (generation === loadGenerationRef.current) {
              component.autoView(
                showNeighborhood && loadedContextSelection ? loadedContextSelection : loadedResidueSelection,
                350
              );
            }
          }, 550);
        }
      })
      .catch((error: unknown) => {
        if (generation === loadGenerationRef.current) {
          setStatus(errorMessage(error, 'Could not load the curated structure'));
        }
      });

    return () => {
      ++loadGenerationRef.current;
    };
    // modelKey captures only properties that require a different coordinate model or selection.
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [modelKey, stageReady]);

  useEffect(() => {
    const component = componentRef.current;
    if (!component) {
      return;
    }
    component.removeAllRepresentations();
    addProteinRepresentations(
      component,
      mode,
      colorMode,
      residueSelection,
      contextSelection,
      showNeighborhood,
      showAssembly,
      chain,
      variant
    );
  }, [chain, colorMode, contextSelection, loadedModel, mode, residueSelection, showAssembly, showNeighborhood, variant]);

  useEffect(() => {
    stageRef.current?.setSpin(spinning);
  }, [spinning]);

  function focusResidue() {
    if (!residueSelection) {
      setStatus('This variant has no curated residue coordinate');
      return;
    }
    componentRef.current?.autoView(showNeighborhood && contextSelection ? contextSelection : residueSelection, 500);
  }

  function resetView() {
    componentRef.current?.autoView('all', 500);
  }

  function exportPng() {
    const stage = stageRef.current;
    if (!stage || !componentRef.current) {
      setStatus('Load a structure before exporting an image');
      return;
    }
    stage
      .makeImage({ factor: 2, antialias: true, trim: true, transparent: false })
      .then((blob) => {
        const url = URL.createObjectURL(blob);
        const link = document.createElement('a');
        link.href = url;
        link.download = `${safeFileName(variant.protein ?? variant.gene ?? 'protein')}-${residueIndex ?? 'model'}.png`;
        link.click();
        window.setTimeout(() => URL.revokeObjectURL(url), 0);
      })
      .catch((error: unknown) => setStatus(errorMessage(error, 'PNG export failed')));
  }

  return (
    <div className="min-w-0 overflow-hidden rounded-md border border-line bg-shell">
      <div className="flex flex-wrap items-center justify-between gap-3 border-b border-line bg-panel2 px-3 py-2">
        <div className="min-w-0 text-sm">
          <div className="truncate font-semibold">{variant.protein ?? variant.gene ?? 'Protein model'}</div>
          <div className="truncate text-xs text-muted" role="status" aria-live="polite">
            {loadedModel?.label ?? structureId ?? 'No structure'} | {status}
          </div>
        </div>
        <div className="flex flex-wrap gap-2" aria-label="Structure representation">
          {MODES.map((item) => (
            <TextButton key={item.value} active={mode === item.value} onClick={() => setMode(item.value)}>
              {item.label}
            </TextButton>
          ))}
        </div>
      </div>
      <div className="relative">
        <div ref={hostRef} className="h-[460px] w-full" aria-label="Interactive protein structure viewer" />
        {!componentRef.current && (
          <div className="pointer-events-none absolute inset-0 grid place-items-center bg-shell/40 px-8 text-center text-sm text-muted">
            {status}
          </div>
        )}
      </div>
      <div className="grid gap-3 border-t border-line bg-panel px-3 py-3 lg:grid-cols-[minmax(0,1fr)_auto]">
        <div className="grid gap-2 text-xs text-muted sm:grid-cols-3">
          <div>
            <span className="font-semibold text-text">Model</span>
            <span className="ml-2">{loadedModel?.coordinateDescription ?? 'not loaded'}</span>
          </div>
          <div>
            <span className="font-semibold text-amber">Orange atoms</span>
            <span className="ml-2">mapped variant residue</span>
          </div>
          <div>
            <span className="font-semibold text-aqua">Context</span>
            <span className="ml-2">5 Å neighborhood around {variant.residue ?? residueSelection ?? 'unmapped residue'}</span>
          </div>
        </div>
        <div className="flex flex-wrap items-center gap-2">
          <select
            value={colorMode}
            onChange={(event) => setColorMode(event.target.value as ColorMode)}
            className="h-8 rounded-md border border-line bg-panel2 px-2 text-xs text-text"
            aria-label="Protein color scheme"
          >
            {COLORS.map((item) => (
              <option key={item.value} value={item.value}>
                {item.value === 'confidence' && loadedModel?.confidenceSemantics === 'bfactor'
                  ? 'Experimental B-factor'
                  : item.label}
              </option>
            ))}
          </select>
          <label className="inline-flex h-8 items-center gap-2 rounded-md border border-line px-2 text-xs text-muted">
            <input
              type="checkbox"
              checked={showNeighborhood}
              onChange={(event) => setShowNeighborhood(event.target.checked)}
            />
            5 Å context
          </label>
          <label className="inline-flex h-8 items-center gap-2 rounded-md border border-line px-2 text-xs text-muted">
            <input
              type="checkbox"
              checked={showAssembly}
              onChange={(event) => setShowAssembly(event.target.checked)}
            />
            Full complex
          </label>
        </div>
      </div>
      <div className="flex flex-wrap items-center justify-between gap-2 border-t border-line bg-panel2 px-3 py-2">
        <div className="min-w-0 text-xs text-muted">
          Selection: <span className="font-mono text-amber">{residueSelection ?? 'not mapped'}</span>
          {loadedModel && (
            <a
              href={loadedModel.entryUrl}
              target="_blank"
              rel="noreferrer"
              className="ml-3 inline-flex items-center gap-1 text-aqua hover:underline"
            >
              Source <ExternalLink className="h-3 w-3" aria-hidden />
            </a>
          )}
        </div>
        <div className="flex flex-wrap gap-2">
          <IconButton label="Focus residue" onClick={focusResidue} icon={<Focus className="h-4 w-4" />} />
          <IconButton label="Reset view" onClick={resetView} icon={<RotateCcw className="h-4 w-4" />} />
          <IconButton
            label={spinning ? 'Pause spin' : 'Spin'}
            onClick={() => setSpinning((current) => !current)}
            icon={spinning ? <Pause className="h-4 w-4" /> : <Play className="h-4 w-4" />}
          />
          <IconButton label="Export PNG" onClick={exportPng} icon={<Download className="h-4 w-4" />} />
        </div>
      </div>
    </div>
  );
}

function TextButton({
  active,
  children,
  onClick
}: {
  active: boolean;
  children: string;
  onClick: () => void;
}) {
  return (
    <button
      type="button"
      onClick={onClick}
      aria-pressed={active}
      className={[
        'rounded-md border px-2.5 py-1.5 text-xs font-semibold',
        active ? 'border-aqua text-aqua' : 'border-line text-muted hover:bg-panel'
      ].join(' ')}
    >
      {children}
    </button>
  );
}

function IconButton({ icon, label, onClick }: { icon: JSX.Element; label: string; onClick: () => void }) {
  return (
    <button
      type="button"
      title={label}
      aria-label={label}
      onClick={onClick}
      className="inline-flex h-8 w-8 items-center justify-center rounded-md border border-line text-muted hover:border-aqua hover:text-aqua"
    >
      {icon}
    </button>
  );
}

function addProteinRepresentations(
  component: StructureComponent,
  mode: RepresentationMode,
  colorMode: ColorMode,
  residueSelection: string | undefined,
  contextSelection: string | undefined,
  showNeighborhood: boolean,
  showAssembly: boolean,
  chain: string,
  variant: ProteinViewerVariant
) {
  const color = colorMode === 'confidence' ? 'bfactor' : colorMode === 'secondary' ? 'sstruc' : '#38bdf8';
  const modelSelection = showAssembly ? 'all' : `:${chain}`;
  if (mode === 'surface') {
    component.addRepresentation('surface', {
      sele: modelSelection,
      color,
      opacity: 0.28,
      side: 'front',
      probeRadius: 1.4,
      surfaceType: 'av'
    } as never);
  } else if (mode === 'cartoon') {
    component.addRepresentation('cartoon', {
      sele: modelSelection,
      color,
      aspectRatio: 4.5,
      radius: 0.62,
      opacity: 0.98
    } as never);
  } else if (mode === 'tube') {
    component.addRepresentation('tube', {
      sele: modelSelection,
      color,
      radius: 0.7,
      opacity: 0.94
    } as never);
  } else {
    component.addRepresentation('backbone', {
      sele: modelSelection,
      color,
      radius: 0.38,
      opacity: 0.96
    } as never);
  }

  if (!residueSelection) {
    return;
  }
  if (showNeighborhood && contextSelection) {
    component.addRepresentation('licorice', {
      sele: `${contextSelection} and not (${residueSelection})`,
      color: 'element',
      radiusScale: 0.65,
      opacity: 0.72
    } as never);
  }
  component.addRepresentation('ball+stick', {
    sele: residueSelection,
    color: 'element',
    multipleBond: true,
    radiusScale: 2.1,
    aspectRatio: 1.4
  } as never);
  component.addRepresentation('spacefill', {
    sele: residueSelection,
    color: '#f97316',
    radiusScale: 0.78,
    opacity: 0.82
  } as never);
  component.addRepresentation('label', {
    sele: `${residueSelection} and .CA`,
    labelType: 'text',
    labelText: [variant.residue ?? variant.label],
    labelGrouping: 'atom',
    color: '#f8fafc',
    fixedSize: true,
    showBackground: true,
    backgroundColor: '#08111f',
    backgroundOpacity: 0.86,
    backgroundMargin: 2,
    yOffset: 1.2,
    zOffset: 2
  } as never);
}

async function loadProteinStructure(
  stage: Stage,
  structureId: string,
  mappedChain: string
): Promise<{ component: StructureComponent } & LoadedModel> {
  const pdbId = normalizePdbId(structureId);
  if (pdbId) {
    return loadFirstAvailable(stage, rcsbCandidates(pdbId, mappedChain));
  }
  throw new Error(`Unsupported structure identifier: ${structureId}. Only curated four-character PDB IDs are enabled.`);
}

async function loadFirstAvailable(
  stage: Stage,
  candidates: StructureSourceCandidate[]
): Promise<{ component: StructureComponent } & LoadedModel> {
  const failures: string[] = [];

  for (const candidate of candidates) {
    try {
      const component = await stage.loadFile(candidate.url, {
        ext: candidate.ext,
        defaultRepresentation: false
      });
      if (!component || !('structure' in component)) {
        failures.push(`${candidate.label}: no structure returned`);
        continue;
      }
      return {
        component: component as StructureComponent,
        label: candidate.label,
        entryUrl: candidate.entryUrl,
        chain: candidate.chain,
        coordinateDescription: candidate.coordinateDescription,
        confidenceSemantics: candidate.confidenceSemantics
      };
    } catch (error: unknown) {
      failures.push(`${candidate.label}: ${errorMessage(error, 'request failed')}`);
    }
  }

  throw new Error(
    `Curated RCSB structure unavailable. ${failures.join('; ')}. ` +
      'Restore the bundled cache with scripts/fetch_protein_structures.sh.'
  );
}

function parseResidueIndex(residue?: string): number | undefined {
  const match = residue?.match(/(\d+)/);
  if (!match) {
    return undefined;
  }
  const value = Number(match[1]);
  return Number.isSafeInteger(value) && value > 0 ? value : undefined;
}

function sanitizeChain(chain: string): string {
  return /^[A-Za-z0-9]$/.test(chain) ? chain : 'A';
}

function safeFileName(value: string): string {
  return value.replace(/[^A-Za-z0-9._-]+/g, '_').replace(/^_+|_+$/g, '') || 'protein';
}

function errorMessage(error: unknown, fallback: string): string {
  return error instanceof Error && error.message ? error.message : fallback;
}
