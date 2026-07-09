use anyhow::{Context, Result};
use serde_json::Value;

const REPORT_RUNTIME: &str = r#"
class MitoCircos {
  constructor(container, data, options = {}) {
    this.container = container;
    this.data = data;
    this.options = Object.assign({ padding: 46 }, options);
    this.layers = { genes: true, coverage: true, clusters: true, svs: true, snps: true };
    this.palette = ['#5eead4', '#38bdf8', '#f472b6', '#facc15', '#a78bfa', '#fb7185'];
  }

  render() {
    this.container.innerHTML = '';
    const width = Math.max(420, this.container.clientWidth || 720);
    const height = Math.min(760, Math.max(500, width * 0.72));
    const radius = Math.min(width, height) / 2 - this.options.padding;
    const cx = width / 2;
    const cy = height / 2;
    const ns = 'http://www.w3.org/2000/svg';
    const svg = document.createElementNS(ns, 'svg');
    svg.setAttribute('viewBox', `0 0 ${width} ${height}`);
    svg.setAttribute('role', 'img');
    svg.setAttribute('aria-label', 'Circular mitochondrial genome visualization');
    this.container.appendChild(svg);
    this.svg = svg;

    const root = document.createElementNS(ns, 'g');
    root.setAttribute('transform', `translate(${cx} ${cy})`);
    svg.appendChild(root);
    this.root = root;

    this.drawTrack(radius + 8, '#263244', 1.5);
    if (this.layers.coverage) this.drawCoverage(radius - 12);
    if (this.layers.genes) this.drawGenes(radius - 48);
    if (this.layers.clusters) this.drawClusters(radius - 90);
    if (this.layers.svs) this.drawSvs(radius - 128);
    if (this.layers.snps) this.drawSnps(radius - 112);
    this.drawCenter();
    return this;
  }

  updateData(data) {
    this.data = data;
    return this.render();
  }

  destroy() {
    this.container.innerHTML = '';
  }

  exportSVG() {
    return new XMLSerializer().serializeToString(this.svg);
  }

  getActiveLayers() {
    return Object.keys(this.layers).filter((key) => this.layers[key]);
  }

  angle(position) {
    const length = this.data.metadata?.reference_length || 16569;
    return (Number(position || 1) - 1) / length * Math.PI * 2 - Math.PI / 2;
  }

  point(position, radius) {
    const angle = this.angle(position);
    return [Math.cos(angle) * radius, Math.sin(angle) * radius];
  }

  arcPath(start, end, radius, sweep = 1) {
    const [x1, y1] = this.point(start, radius);
    const [x2, y2] = this.point(end, radius);
    const delta = Math.abs(this.angle(end) - this.angle(start));
    const large = delta > Math.PI ? 1 : 0;
    return `M ${x1.toFixed(2)} ${y1.toFixed(2)} A ${radius} ${radius} 0 ${large} ${sweep} ${x2.toFixed(2)} ${y2.toFixed(2)}`;
  }

  node(name, attrs = {}) {
    const element = document.createElementNS('http://www.w3.org/2000/svg', name);
    for (const [key, value] of Object.entries(attrs)) {
      element.setAttribute(key, String(value));
    }
    return element;
  }

  drawTrack(radius, color, width) {
    this.root.appendChild(this.node('circle', {
      r: radius,
      fill: 'none',
      stroke: color,
      'stroke-width': width
    }));
  }

  drawCoverage(radius) {
    const coverage = this.data.coverage || [];
    const maxDepth = Math.max(1, ...coverage.map((bin) => Number(bin.depth || 0)));
    for (const bin of coverage) {
      const depth = Number(bin.depth || 0);
      const width = 2 + depth / maxDepth * 34;
      const path = this.node('path', {
        d: this.arcPath(bin.start, bin.end, radius),
        fill: 'none',
        stroke: '#67e8f9',
        'stroke-opacity': 0.26 + depth / maxDepth * 0.60,
        'stroke-width': width,
        'stroke-linecap': 'butt'
      });
      this.root.appendChild(path);
    }
  }

  drawGenes(radius) {
    const genes = this.data.genes || [];
    genes.forEach((gene, index) => {
      const color = this.palette[index % this.palette.length];
      const path = this.node('path', {
        d: this.arcPath(gene.start, gene.end, radius),
        fill: 'none',
        stroke: color,
        'stroke-width': 18,
        'stroke-linecap': 'round',
        tabindex: 0
      });
      path.addEventListener('mouseenter', () => this.showTooltip(`${gene.name} ${gene.start}-${gene.end}`));
      path.addEventListener('mouseleave', () => this.showTooltip(''));
      path.addEventListener('click', () => {
        document.querySelector('#selection').textContent = `${gene.name} (${gene.biotype})`;
      });
      this.root.appendChild(path);

      if (index % 2 === 0) {
        const [x, y] = this.point((gene.start + gene.end) / 2, radius + 24);
        const label = this.node('text', {
          x,
          y,
          fill: '#cbd5e1',
          'font-size': 10,
          'text-anchor': 'middle',
          'dominant-baseline': 'middle'
        });
        label.textContent = gene.name.replace('MT-', '');
        this.root.appendChild(label);
      }
    });
  }

  drawClusters(radius) {
    const clusters = this.data.clusters || [];
    const maxSize = Math.max(1, ...clusters.map((cluster) => Number(cluster.size || 0)));
    clusters.forEach((cluster, index) => {
      const offset = index * 13;
      const path = this.node('circle', {
        r: radius - offset,
        fill: 'none',
        stroke: this.palette[index % this.palette.length],
        'stroke-width': 3 + Number(cluster.size || 0) / maxSize * 8,
        'stroke-opacity': 0.78
      });
      path.addEventListener('click', () => {
        document.querySelector('#selection').textContent = `${cluster.label}: ${cluster.size} reads`;
        if (window.renderClusterDetails) window.renderClusterDetails(cluster.id);
      });
      this.root.appendChild(path);
    });
  }

  drawSvs(radius) {
    for (const sv of this.data.svs || []) {
      const [x1, y1] = this.point(sv.start, radius);
      const [x2, y2] = this.point(sv.end, radius);
      const stroke = sv.known_event ? '#f97316' : '#e879f9';
      const width = Math.min(10, 2 + (sv.supporting_reads?.length || 1));
      const path = this.node('path', {
        d: `M ${x1.toFixed(1)} ${y1.toFixed(1)} Q 0 0 ${x2.toFixed(1)} ${y2.toFixed(1)}`,
        fill: 'none',
        stroke,
        'stroke-width': width,
        'stroke-opacity': 0.66
      });
      path.addEventListener('mouseenter', () => {
        path.setAttribute('stroke', '#f8fafc');
        path.setAttribute('stroke-width', String(Math.min(16, width + 5)));
        path.setAttribute('stroke-opacity', '1');
        this.showTooltip(`${sv.type} ${sv.start}-${sv.end}`);
      });
      path.addEventListener('mouseleave', () => {
        path.setAttribute('stroke', stroke);
        path.setAttribute('stroke-width', String(width));
        path.setAttribute('stroke-opacity', '0.66');
        this.showTooltip('');
      });
      path.addEventListener('click', () => {
        document.querySelector('#selection').textContent = `${sv.type} ${sv.start}-${sv.end}, n=${sv.supporting_reads?.length || 0}`;
      });
      this.root.appendChild(path);
    }
  }

  drawSnps(radius) {
    const seen = new Map();
    for (const read of this.data.reads || []) {
      if (read.filtered_numt) continue;
      for (const snp of read.snps || []) {
        const key = `${snp.position}:${snp.alt}`;
        seen.set(key, snp);
      }
    }
    for (const snp of seen.values()) {
      const [x, y] = this.point(snp.position, radius);
      this.root.appendChild(this.node('circle', {
        cx: x,
        cy: y,
        r: 3.2,
        fill: snp.alt === 'A' ? '#22c55e' : snp.alt === 'C' ? '#06b6d4' : snp.alt === 'G' ? '#f59e0b' : '#ef4444',
        stroke: '#0f172a',
        'stroke-width': 1
      }));
    }
  }

  drawCenter() {
    const passed = this.data.filter_stats?.passed_reads ?? 0;
    const filtered = this.data.filter_stats?.numt_filtered_reads ?? 0;
    const title = this.node('text', {
      x: 0,
      y: -8,
      fill: '#f8fafc',
      'font-size': 24,
      'font-weight': 700,
      'text-anchor': 'middle'
    });
    title.textContent = this.data.metadata?.sample || 'mtDNA';
    this.root.appendChild(title);
    const subtitle = this.node('text', {
      x: 0,
      y: 20,
      fill: '#94a3b8',
      'font-size': 13,
      'text-anchor': 'middle'
    });
    subtitle.textContent = `${passed} reads passed, ${filtered} NUMT-filtered`;
    this.root.appendChild(subtitle);
  }

  showTooltip(text) {
    const tooltip = document.querySelector('#tooltip');
    if (tooltip) tooltip.textContent = text;
  }
}
window.MitoCircos = MitoCircos;
"#;

pub fn render_report(json: &str) -> Result<String> {
    let data: Value = serde_json::from_str(json).context("analysis returned invalid JSON")?;
    let data_literal = serde_json::to_string(&data)?
        .replace("</", "<\\/")
        .replace("<!--", "<\\!--");
    let pretty_json = serde_json::to_string_pretty(&data)?;

    Ok(include_str!("../templates/report.html")
        .replace("{{DATA}}", &data_literal)
        .replace("{{RUNTIME_JS}}", REPORT_RUNTIME)
        .replace("{{PRETTY_JSON}}", &escape_html(&pretty_json)))
}

fn escape_html(value: &str) -> String {
    value
        .replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
        .replace('\'', "&#39;")
}
