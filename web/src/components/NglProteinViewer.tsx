import { Download, Focus, Pause, Play, RotateCcw } from 'lucide-react';
import { useEffect, useMemo, useRef, useState } from 'react';
import type { ProteinStructureMapping } from '@mito-architect/visualization-lib';
import type { Component, Stage, StructureComponent } from 'ngl';

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

interface NglProteinViewerProps {
  variant: ProteinViewerVariant;
}

const MODES: Array<{ value: RepresentationMode; label: string }> = [
  { value: 'cartoon', label: 'Cartoon' },
  { value: 'tube', label: 'Tube' },
  { value: 'backbone', label: 'Backbone' },
  { value: 'surface', label: 'Surface' }
];

export default function NglProteinViewer({ variant }: NglProteinViewerProps) {
  const hostRef = useRef<HTMLDivElement>(null);
  const stageRef = useRef<Stage>();
  const componentRef = useRef<StructureComponent>();
  const [mode, setMode] = useState<RepresentationMode>('cartoon');
  const [spinning, setSpinning] = useState(false);
  const [status, setStatus] = useState('Loading structure');
  const [stageReady, setStageReady] = useState(false);
  const residueIndex = variant.structure?.residue_index ?? parseResidueIndex(variant.residue) ?? 1;
  const chain = variant.structure?.chain ?? 'A';
  const residueSelection = `${residueIndex}:${chain}`;
  const fullSelection = 'all';
  const pdb = useMemo(() => buildLocalProteinModel(variant), [variant]);

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
      .catch(() => setStatus('NGL failed to load'));

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

  useEffect(() => {
    const stage = stageRef.current;
    if (!stage || !stageReady) {
      return;
    }

    let cancelled = false;
    setStatus('Loading structure');
    stage.setSpin(false);
    stage.removeAllComponents();

    loadProteinStructure(stage, variant, pdb)
      .then(({ component, source }) => {
        if (!component) {
          setStatus('No structure component returned');
          return;
        }
        if (cancelled) {
          stage.removeComponent(component);
          return;
        }

        const structureComponent = component as StructureComponent;
        componentRef.current = structureComponent;
        addProteinRepresentations(structureComponent, mode, fullSelection, residueSelection, variant);
        structureComponent.autoView(fullSelection, 500);
        window.setTimeout(() => structureComponent.autoView(residueSelection, 350), 450);
        stage.setSpin(spinning);
        setStatus(source === 'alphafold' ? 'AlphaFold model' : 'Local fallback model');
      })
      .catch((error: unknown) => {
        setStatus(error instanceof Error ? error.message : 'Could not load structure');
      });

    return () => {
      cancelled = true;
    };
  }, [fullSelection, mode, pdb, residueSelection, spinning, stageReady, variant]);

  function focusResidue() {
    componentRef.current?.autoView(residueSelection, 500);
  }

  function resetView() {
    componentRef.current?.autoView(fullSelection, 500);
  }

  function toggleSpin() {
    setSpinning((current) => {
      const next = !current;
      stageRef.current?.setSpin(next);
      return next;
    });
  }

  function exportPng() {
    stageRef.current
      ?.makeImage({ factor: 2, antialias: true, trim: true })
      .then((blob) => {
        const url = URL.createObjectURL(blob);
        const link = document.createElement('a');
        link.href = url;
        link.download = `${variant.protein ?? variant.gene ?? 'protein'}-${residueIndex}.png`;
        link.click();
        URL.revokeObjectURL(url);
      })
      .catch(() => setStatus('PNG export failed'));
  }

  return (
    <div className="min-w-0 overflow-hidden rounded-md border border-line bg-shell">
      <div className="flex flex-wrap items-center justify-between gap-2 border-b border-line bg-panel2 px-3 py-2">
        <div className="min-w-0 text-sm">
          <div className="truncate font-semibold">{variant.protein ?? variant.gene ?? 'Protein model'}</div>
          <div className="truncate text-xs text-muted">
            {variant.structure?.structure_id ?? 'local PDB model'} | {status}
          </div>
        </div>
        <div className="flex flex-wrap gap-2">
          {MODES.map((item) => (
            <button
              key={item.value}
              type="button"
              onClick={() => setMode(item.value)}
              className={[
                'rounded-md border px-2.5 py-1.5 text-xs font-semibold',
                mode === item.value ? 'border-aqua text-aqua' : 'border-line text-muted hover:bg-panel'
              ].join(' ')}
            >
              {item.label}
            </button>
          ))}
        </div>
      </div>
      <div ref={hostRef} className="h-[420px] w-full" />
      <div className="border-t border-line bg-panel px-3 py-2">
        <div className="grid gap-2 text-xs text-muted sm:grid-cols-3">
          <div>
            <span className="font-semibold text-text">Whole model</span>
            <span className="ml-2">shown by the selected representation</span>
          </div>
          <div>
            <span className="font-semibold text-amber">Orange atoms</span>
            <span className="ml-2">mutated residue position</span>
          </div>
          <div>
            <span className="font-semibold text-aqua">Residue</span>
            <span className="ml-2">
              {variant.residue ?? residueSelection} in {variant.protein ?? variant.gene}
            </span>
          </div>
        </div>
      </div>
      <div className="flex flex-wrap items-center justify-between gap-2 border-t border-line bg-panel2 px-3 py-2">
        <div className="min-w-0 text-xs text-muted">
          Highlight: <span className="font-semibold text-amber">{variant.residue ?? residueSelection}</span>
          <span className="mx-2">|</span>
          Selection: <span className="font-mono">{residueSelection}</span>
        </div>
        <div className="flex flex-wrap gap-2">
          <IconButton label="Focus residue" onClick={focusResidue} icon={<Focus className="h-4 w-4" />} />
          <IconButton label="Reset view" onClick={resetView} icon={<RotateCcw className="h-4 w-4" />} />
          <IconButton
            label={spinning ? 'Pause spin' : 'Spin'}
            onClick={toggleSpin}
            icon={spinning ? <Pause className="h-4 w-4" /> : <Play className="h-4 w-4" />}
          />
          <IconButton label="Export PNG" onClick={exportPng} icon={<Download className="h-4 w-4" />} />
        </div>
      </div>
    </div>
  );
}

function IconButton({
  icon,
  label,
  onClick
}: {
  icon: JSX.Element;
  label: string;
  onClick: () => void;
}) {
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
  fullSelection: string,
  residueSelection: string,
  variant: ProteinViewerVariant
) {
  if (mode === 'surface') {
    component.addRepresentation('surface', {
      sele: fullSelection,
      color: '#38bdf8',
      opacity: 0.32,
      side: 'front'
    } as never);
    component.addRepresentation('licorice', {
      sele: fullSelection,
      color: 'element',
      radius: 0.12,
      opacity: 0.55
    } as never);
  } else if (mode === 'cartoon') {
    component.addRepresentation('cartoon', {
      sele: fullSelection,
      color: 'sstruc',
      aspectRatio: 4.5,
      radius: 0.62,
      opacity: 0.98
    } as never);
    component.addRepresentation('trace', {
      sele: fullSelection,
      color: '#38bdf8',
      radius: 0.18,
      opacity: 0.5
    } as never);
  } else if (mode === 'tube') {
    component.addRepresentation('tube', {
      sele: fullSelection,
      color: '#5eead4',
      radius: 0.7,
      opacity: 0.94
    } as never);
  } else {
    component.addRepresentation('backbone', {
      sele: fullSelection,
      color: '#a78bfa',
      radius: 0.38,
      opacity: 0.96
    } as never);
    component.addRepresentation('line', {
      sele: fullSelection,
      color: 'element',
      opacity: 0.45
    } as never);
  }

  component.addRepresentation('ball+stick', {
    sele: residueSelection,
    color: 'element',
    multipleBond: true,
    radiusScale: 2.2,
    aspectRatio: 1.4
  } as never);
  component.addRepresentation('spacefill', {
    sele: residueSelection,
    color: '#f97316',
    radiusScale: 0.82,
    opacity: 0.88
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
    backgroundOpacity: 0.82,
    backgroundMargin: 2,
    yOffset: 1.2,
    zOffset: 2
  } as never);
}

async function loadProteinStructure(
  stage: Stage,
  variant: ProteinViewerVariant,
  fallbackPdb: string
): Promise<{ component: Component | void; source: 'alphafold' | 'fallback' }> {
  const alphafoldUrl = alphaFoldPdbUrl(variant.structure?.structure_id);
  if (alphafoldUrl) {
    try {
      const component = await stage.loadFile(alphafoldUrl, {
        ext: 'pdb',
        defaultRepresentation: false
      });
      return { component, source: 'alphafold' };
    } catch {
      // Offline reports and restricted networks still get an inspectable model.
    }
  }

  const blob = new Blob([fallbackPdb], { type: 'text/plain' });
  const component = await stage.loadFile(blob, { ext: 'pdb', defaultRepresentation: false });
  return { component, source: 'fallback' };
}

function alphaFoldPdbUrl(structureId?: string): string | undefined {
  if (!structureId?.startsWith('AF-')) {
    return undefined;
  }
  const normalized = structureId.endsWith('-model_v4') ? structureId : `${structureId}-model_v4`;
  return `https://alphafold.ebi.ac.uk/files/${normalized}.pdb`;
}

function buildLocalProteinModel(variant: ProteinViewerVariant): string {
  const residueIndex = variant.structure?.residue_index ?? parseResidueIndex(variant.residue) ?? 156;
  const chain = (variant.structure?.chain ?? 'A').slice(0, 1);
  const start = Math.max(1, residueIndex - 38);
  const end = residueIndex + 38;
  const mutant = parseMutantResidue(variant.residue);
  const lines = [
    `HEADER    MITO-ARCHITECT LOCAL PROTEIN MODEL`,
    `TITLE     ${variant.protein ?? variant.gene ?? 'MITOCHONDRIAL PROTEIN'} ${variant.residue ?? variant.label}`,
    helixLine(1, 'H1', chain, start, end)
  ];
  let serial = 1;

  for (let resno = start; resno <= end; ++resno) {
    const offset = resno - start;
    const angle = offset * 1.7453292519943295;
    const bend = Math.sin(offset / 9) * 4.2;
    const x = Math.cos(angle) * 3.1 + bend;
    const y = Math.sin(angle) * 3.1;
    const z = offset * 1.48 - 55;
    const resName = resno === residueIndex ? mutant : residueNameForOffset(offset);

    lines.push(atomLine(serial++, 'N', resName, chain, resno, x - 1.25, y + 0.45, z - 0.52, 'N'));
    lines.push(atomLine(serial++, 'CA', resName, chain, resno, x, y, z, 'C'));
    lines.push(atomLine(serial++, 'C', resName, chain, resno, x + 1.23, y - 0.35, z + 0.52, 'C'));
    lines.push(atomLine(serial++, 'O', resName, chain, resno, x + 1.92, y - 1.08, z + 0.24, 'O'));
    lines.push(atomLine(serial++, 'CB', resName, chain, resno, x - 0.28, y + 1.46, z + 0.86, 'C'));
  }

  lines.push('TER');
  lines.push('END');
  return `${lines.join('\n')}\n`;
}

function atomLine(
  serial: number,
  atom: string,
  resName: string,
  chain: string,
  resno: number,
  x: number,
  y: number,
  z: number,
  element: string
): string {
  return [
    'ATOM  ',
    serial.toString().padStart(5, ' '),
    ' ',
    atom.padEnd(4, ' '),
    ' ',
    resName.padStart(3, ' '),
    ' ',
    chain,
    resno.toString().padStart(4, ' '),
    '    ',
    x.toFixed(3).padStart(8, ' '),
    y.toFixed(3).padStart(8, ' '),
    z.toFixed(3).padStart(8, ' '),
    '  1.00',
    ' 20.00',
    '           ',
    element.padStart(2, ' ')
  ].join('');
}

function helixLine(serial: number, id: string, chain: string, start: number, end: number): string {
  return `HELIX  ${serial.toString().padStart(3, ' ')} ${id.padEnd(3, ' ')} ALA ${chain}${start
    .toString()
    .padStart(4, ' ')}  ALA ${chain}${end.toString().padStart(4, ' ')}  1${(end - start + 1)
    .toString()
    .padStart(36, ' ')}`;
}

function parseResidueIndex(residue?: string): number | undefined {
  const match = residue?.match(/(\d+)/);
  return match ? Number(match[1]) : undefined;
}

function parseMutantResidue(residue?: string): string {
  const match = residue?.match(/[A-Z][a-z]{2}\d+([A-Z][a-z]{2})$/);
  return match ? match[1].toUpperCase() : 'ARG';
}

function residueNameForOffset(offset: number): string {
  const residues = ['ALA', 'LEU', 'ILE', 'VAL', 'GLY', 'SER', 'THR', 'PHE', 'TYR', 'MET'];
  return residues[offset % residues.length];
}
