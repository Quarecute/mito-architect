import react from '@vitejs/plugin-react';
import { defineConfig } from 'vite';

export default defineConfig({
  plugins: [react()],
  resolve: {
    alias: {
      '@mito-architect/visualization-lib': new URL(
        '../visualization-lib/src/index.ts',
        import.meta.url
      ).pathname
    }
  },
  server: {
    host: '127.0.0.1',
    port: 5173
  },
  build: {
    chunkSizeWarningLimit: 1400
  }
});
