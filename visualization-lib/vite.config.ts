import { defineConfig } from 'vite';

export default defineConfig({
  build: {
    lib: {
      entry: 'src/index.ts',
      name: 'MitoCircos',
      fileName: (format) => (format === 'umd' ? 'mito-circos.umd.cjs' : 'mito-circos.js'),
      formats: ['es', 'umd']
    },
    rollupOptions: {
      external: [],
      output: {
        globals: {
          d3: 'd3'
        }
      }
    },
    sourcemap: true
  }
});

