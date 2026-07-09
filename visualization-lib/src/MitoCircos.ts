import * as d3 from 'd3';
import type {
  ClusterSummary,
  GeneAnnotation,
  LayerName,
  LayerState,
  MitoAnalysisData,
  MitoCircosOptions,
  ReadFeature,
  SnpCall,
  StructuralVariant
} from './types';

const DEFAULT_LAYERS: LayerState = {
  genes: true,
  coverage: true,
  clusters: true,
  svs: true,
  snps: true
};

const DEFAULT_PALETTE = ['#5eead4', '#38bdf8', '#f472b6', '#facc15', '#a78bfa', '#fb7185'];

export class MitoCircos {
  private container: HTMLElement;
  private data: MitoAnalysisData;
  private options: MitoCircosOptions;
  private layers: LayerState;
  private svg?: d3.Selection<SVGSVGElement, unknown, null, undefined>;
  private plot?: d3.Selection<SVGGElement, unknown, null, undefined>;
  private tooltip?: HTMLDivElement;
  private width = 760;
  private height = 620;
  private radius = 250;
  private angleScale!: d3.ScaleLinear<number, number>;

  constructor(container: HTMLElement, data: MitoAnalysisData, options: MitoCircosOptions = {}) {
    this.container = container;
    this.data = data;
    this.options = options;
    this.layers = { ...DEFAULT_LAYERS, ...options.activeLayers };
  }

  /** Render the circular mtDNA plot into the configured container. */
  render(): this {
    this.destroy();
    this.measure();
    this.angleScale = d3
      .scaleLinear()
      .domain([1, Math.max(2, this.data.metadata.reference_length + 1)])
      .range([0, Math.PI * 2]);

    const root = d3
      .select(this.container)
      .append('div')
      .attr('class', `mito-circos mito-circos-${this.options.theme ?? 'dark'}`);

    if (this.options.showControls ?? true) {
      this.renderControls(root);
    }

    const frame = root.append('div').attr('class', 'mito-circos-frame');
    this.svg = frame
      .append('svg')
      .attr('viewBox', `0 0 ${this.width} ${this.height}`)
      .attr('role', 'img')
      .attr('aria-label', `Circular mtDNA plot for ${this.data.metadata.sample}`)
      .attr('class', 'mito-circos-svg');

    this.tooltip = frame.append('div').attr('class', 'mito-circos-tooltip').node() ?? undefined;

    this.plot = this.svg
      .append('g')
      .attr('transform', `translate(${this.width / 2}, ${this.height / 2})`);

    this.svg.call(
      d3
        .zoom<SVGSVGElement, unknown>()
        .scaleExtent([0.7, 6])
        .on('zoom', (event) => {
          this.plot?.attr(
            'transform',
            `translate(${this.width / 2}, ${this.height / 2}) ${event.transform.toString()}`
          );
        })
    );

    this.drawBackbone();
    if (this.layers.coverage) this.drawCoverage();
    if (this.layers.genes) this.drawGenes();
    if (this.layers.clusters) this.drawClusters();
    if (this.layers.svs) this.drawStructuralVariants();
    if (this.layers.snps) this.drawSnps();
    this.drawCenterLabel();
    this.injectStyles();
    return this;
  }

  /** Replace the result payload and redraw all active layers. */
  updateData(data: MitoAnalysisData): this {
    this.data = data;
    return this.render();
  }

  /** Remove generated DOM nodes and release references held by this instance. */
  destroy(): void {
    this.container.innerHTML = '';
    this.svg = undefined;
    this.plot = undefined;
    this.tooltip = undefined;
  }

  /** Serialize the current SVG plot for report export or download buttons. */
  exportSVG(): string {
    if (!this.svg) {
      return '';
    }
    const node = this.svg.node();
    if (!node) {
      return '';
    }
    const clone = node.cloneNode(true) as SVGSVGElement;
    clone.setAttribute('xmlns', 'http://www.w3.org/2000/svg');
    return new XMLSerializer().serializeToString(clone);
  }

  /** Return the names of layers currently visible in the plot. */
  getActiveLayers(): LayerName[] {
    return (Object.keys(this.layers) as LayerName[]).filter((layer) => this.layers[layer]);
  }

  /** Toggle one visualization layer and redraw the plot. */
  setLayer(layer: LayerName, active: boolean): this {
    this.layers[layer] = active;
    return this.render();
  }

  private measure(): void {
    const bounds = this.container.getBoundingClientRect();
    this.width = Math.max(420, this.options.width ?? bounds.width ?? 760);
    this.height = Math.max(460, this.options.height ?? Math.min(780, this.width * 0.72));
    this.radius = Math.min(this.width, this.height) / 2 - 58;
  }

  private renderControls(root: d3.Selection<HTMLDivElement, unknown, null, undefined>): void {
    const labels: Array<[LayerName, string]> = [
      ['genes', 'Genes'],
      ['coverage', 'Coverage'],
      ['clusters', 'Clusters'],
      ['svs', 'SVs'],
      ['snps', 'SNPs']
    ];
    const controls = root.append('div').attr('class', 'mito-circos-controls');
    for (const [layer, label] of labels) {
      const item = controls.append('label').attr('class', 'mito-circos-toggle');
      item
        .append('input')
        .attr('type', 'checkbox')
        .property('checked', this.layers[layer])
        .on('change', (event) => {
          this.layers[layer] = Boolean((event.currentTarget as HTMLInputElement).checked);
          this.render();
        });
      item.append('span').text(label);
    }
  }

  private drawBackbone(): void {
    this.plot
      ?.append('circle')
      .attr('r', this.radius + 18)
      .attr('fill', 'none')
      .attr('stroke', 'var(--mito-border)')
      .attr('stroke-width', 1.5);
  }

  private drawCoverage(): void {
    const maxDepth = Math.max(1, d3.max(this.data.coverage, (d) => d.depth) ?? 1);
    const track = this.plot?.append('g').attr('class', 'mito-circos-coverage');
    const arc = d3
      .arc<unknown>()
      .innerRadius(this.radius - 34)
      .outerRadius((d: unknown) => {
        const bin = d as { depth: number };
        return this.radius - 34 + (bin.depth / maxDepth) * 42;
      })
      .startAngle((d: unknown) => this.angle((d as { start: number }).start))
      .endAngle((d: unknown) => this.angle((d as { end: number }).end + 1));

    track
      ?.selectAll('path')
      .data(this.data.coverage)
      .join('path')
      .attr('d', (d) => arc(d) ?? '')
      .attr('fill', '#67e8f9')
      .attr('fill-opacity', 0.22)
      .attr('stroke', '#0891b2')
      .attr('stroke-opacity', 0.28)
      .on('pointermove', (event, d) =>
        this.showTooltip(event, `Coverage ${d.start}-${d.end}: ${d.depth}`)
      )
      .on('pointerleave', () => this.hideTooltip());
  }

  private drawGenes(): void {
    const track = this.plot?.append('g').attr('class', 'mito-circos-genes');
    const arc = d3
      .arc<GeneAnnotation>()
      .innerRadius(this.radius - 78)
      .outerRadius(this.radius - 54)
      .startAngle((d) => this.angle(d.start))
      .endAngle((d) => this.angle(d.end + 1))
      .padAngle(0.004)
      .cornerRadius(6);

    track
      ?.selectAll('path')
      .data(this.data.genes)
      .join('path')
      .attr('d', (d) => arc(d) ?? '')
      .attr('fill', (_, index) => this.color(index))
      .attr('stroke', 'rgba(8, 17, 31, 0.7)')
      .attr('stroke-width', 1)
      .attr('tabindex', 0)
      .on('pointermove', (event, d) =>
        this.showTooltip(event, `${d.name} ${d.start}-${d.end} ${d.strand}`)
      )
      .on('pointerleave', () => this.hideTooltip())
      .on('click', (_, d) => this.dispatch('mito:gene-select', d));

    const labels = track
      ?.selectAll('text')
      .data(this.data.genes.filter((_, index) => index % 2 === 0))
      .join('text')
      .attr('class', 'mito-circos-label')
      .attr('x', (d) => this.point((d.start + d.end) / 2, this.radius - 38)[0])
      .attr('y', (d) => this.point((d.start + d.end) / 2, this.radius - 38)[1])
      .attr('text-anchor', 'middle')
      .attr('dominant-baseline', 'middle')
      .text((d) => d.name.replace('MT-', ''));

    labels?.attr('transform', (d) => {
      const [x, y] = this.point((d.start + d.end) / 2, this.radius - 38);
      return `rotate(${this.degrees((d.start + d.end) / 2)}, ${x}, ${y})`;
    });
  }

  private drawClusters(): void {
    const maxSize = Math.max(1, d3.max(this.data.clusters, (d) => d.size) ?? 1);
    const track = this.plot?.append('g').attr('class', 'mito-circos-clusters');
    track
      ?.selectAll('circle')
      .data(this.data.clusters)
      .join('circle')
      .attr('r', (_, index) => Math.max(8, this.radius - 110 - index * 14))
      .attr('fill', 'none')
      .attr('stroke', (_, index) => this.color(index))
      .attr('stroke-width', (d) => 2 + (d.size / maxSize) * 8)
      .attr('stroke-opacity', 0.78)
      .on('pointermove', (event, d) => this.showTooltip(event, `${d.label}: ${d.size} reads`))
      .on('pointerleave', () => this.hideTooltip())
      .on('click', (_, d) => this.dispatch('mito:cluster-select', d));
  }

  private drawStructuralVariants(): void {
    const track = this.plot?.append('g').attr('class', 'mito-circos-svs');
    track
      ?.selectAll('path')
      .data(this.data.svs)
      .join('path')
      .attr('d', (d) => this.chordPath(d))
      .attr('fill', 'none')
      .attr('stroke', (d) => (d.known_event ? '#f97316' : '#e879f9'))
      .attr('stroke-width', (d) => Math.min(12, 2 + d.supporting_reads.length))
      .attr('stroke-opacity', 0.68)
      .on('pointerenter', (event, d) => {
        d3.select(event.currentTarget as SVGPathElement)
          .attr('stroke', '#f8fafc')
          .attr('stroke-width', Math.min(16, 5 + d.supporting_reads.length))
          .attr('stroke-opacity', 1);
      })
      .on('pointermove', (event, d) =>
        this.showTooltip(event, `${d.type} ${d.start}-${d.end}, n=${d.supporting_reads.length}`)
      )
      .on('pointerleave', (event, d) => {
        d3.select(event.currentTarget as SVGPathElement)
          .attr('stroke', d.known_event ? '#f97316' : '#e879f9')
          .attr('stroke-width', Math.min(12, 2 + d.supporting_reads.length))
          .attr('stroke-opacity', 0.68);
        this.hideTooltip();
      })
      .on('click', (_, d) => this.dispatch('mito:sv-select', d));
  }

  private drawSnps(): void {
    const snps = this.uniqueSnps(this.data.reads);
    const track = this.plot?.append('g').attr('class', 'mito-circos-snps');
    track
      ?.selectAll('circle')
      .data(snps)
      .join('circle')
      .attr('cx', (d) => this.point(d.position, this.radius - 118)[0])
      .attr('cy', (d) => this.point(d.position, this.radius - 118)[1])
      .attr('r', 3.6)
      .attr('fill', (d) => this.snpColor(d.alt))
      .attr('stroke', 'var(--mito-bg)')
      .attr('stroke-width', 1)
      .on('pointermove', (event, d) => this.showTooltip(event, `${d.position} ${d.ref}>${d.alt}`))
      .on('pointerleave', () => this.hideTooltip())
      .on('click', (_, d) => this.dispatch('mito:snp-select', d));
  }

  private drawCenterLabel(): void {
    const group = this.plot?.append('g').attr('class', 'mito-circos-center');
    group
      ?.append('text')
      .attr('text-anchor', 'middle')
      .attr('y', -12)
      .attr('class', 'mito-circos-title')
      .text(this.data.metadata.sample);
    group
      ?.append('text')
      .attr('text-anchor', 'middle')
      .attr('y', 16)
      .attr('class', 'mito-circos-subtitle')
      .text(`${this.data.filter_stats.passed_reads} reads, ${this.data.clusters.length} clusters`);
  }

  private angle(position: number): number {
    return this.angleScale(Math.max(1, Math.min(this.data.metadata.reference_length, position)));
  }

  private point(position: number, radius: number): [number, number] {
    const angle = this.angle(position) - Math.PI / 2;
    return [Math.cos(angle) * radius, Math.sin(angle) * radius];
  }

  private degrees(position: number): number {
    return (this.angle(position) * 180) / Math.PI;
  }

  private chordPath(sv: StructuralVariant): string {
    const [x1, y1] = this.point(sv.start, this.radius - 132);
    const [x2, y2] = this.point(sv.end, this.radius - 132);
    return `M ${x1.toFixed(2)} ${y1.toFixed(2)} Q 0 0 ${x2.toFixed(2)} ${y2.toFixed(2)}`;
  }

  private uniqueSnps(reads: ReadFeature[]): SnpCall[] {
    const byKey = new Map<string, SnpCall>();
    for (const read of reads) {
      if (read.filtered_numt) continue;
      for (const snp of read.snps) {
        byKey.set(`${snp.position}:${snp.ref}:${snp.alt}`, snp);
      }
    }
    return [...byKey.values()];
  }

  private color(index: number): string {
    const palette = this.options.palette?.length ? this.options.palette : DEFAULT_PALETTE;
    return palette[index % palette.length];
  }

  private snpColor(base: string): string {
    switch (base) {
      case 'A':
        return '#22c55e';
      case 'C':
        return '#06b6d4';
      case 'G':
        return '#f59e0b';
      case 'T':
        return '#ef4444';
      default:
        return '#94a3b8';
    }
  }

  private showTooltip(event: PointerEvent, text: string): void {
    if (!this.tooltip) return;
    const host = this.container.getBoundingClientRect();
    this.tooltip.textContent = text;
    this.tooltip.style.opacity = '1';
    this.tooltip.style.left = `${event.clientX - host.left + 14}px`;
    this.tooltip.style.top = `${event.clientY - host.top + 14}px`;
  }

  private hideTooltip(): void {
    if (!this.tooltip) return;
    this.tooltip.style.opacity = '0';
  }

  private dispatch(name: string, detail: GeneAnnotation | ClusterSummary | StructuralVariant | SnpCall): void {
    this.container.dispatchEvent(new CustomEvent(name, { detail, bubbles: true }));
  }

  private injectStyles(): void {
    if (document.getElementById('mito-circos-styles')) {
      return;
    }
    const style = document.createElement('style');
    style.id = 'mito-circos-styles';
    style.textContent = `
.mito-circos {
  --mito-bg: #08111f;
  --mito-panel: #101a2b;
  --mito-border: #263244;
  --mito-text: #e5edf8;
  --mito-muted: #94a3b8;
  width: 100%;
  color: var(--mito-text);
}
.mito-circos-light {
  --mito-bg: #f8fafc;
  --mito-panel: #ffffff;
  --mito-border: #d7e0eb;
  --mito-text: #172033;
  --mito-muted: #62748e;
}
.mito-circos-frame {
  position: relative;
  overflow: hidden;
  border: 1px solid var(--mito-border);
  border-radius: 8px;
  background: var(--mito-panel);
}
.mito-circos-svg {
  display: block;
  width: 100%;
  height: auto;
  min-height: 420px;
  cursor: grab;
}
.mito-circos-svg:active {
  cursor: grabbing;
}
.mito-circos-controls {
  display: flex;
  flex-wrap: wrap;
  gap: 10px;
  margin-bottom: 10px;
}
.mito-circos-toggle {
  display: inline-flex;
  align-items: center;
  gap: 8px;
  border: 1px solid var(--mito-border);
  border-radius: 8px;
  padding: 7px 10px;
  background: var(--mito-panel);
  font: 600 13px/1.2 ui-sans-serif, system-ui, sans-serif;
}
.mito-circos-label {
  fill: var(--mito-muted);
  font: 600 10px/1 ui-sans-serif, system-ui, sans-serif;
}
.mito-circos-title {
  fill: var(--mito-text);
  font: 750 24px/1 ui-sans-serif, system-ui, sans-serif;
}
.mito-circos-subtitle {
  fill: var(--mito-muted);
  font: 600 13px/1 ui-sans-serif, system-ui, sans-serif;
}
.mito-circos-tooltip {
  position: absolute;
  pointer-events: none;
  opacity: 0;
  z-index: 4;
  max-width: 280px;
  border: 1px solid var(--mito-border);
  border-radius: 8px;
  background: rgba(8, 17, 31, 0.94);
  color: #f8fafc;
  padding: 7px 9px;
  font: 600 12px/1.35 ui-sans-serif, system-ui, sans-serif;
  transition: opacity 120ms ease;
}`;
    document.head.appendChild(style);
  }
}
