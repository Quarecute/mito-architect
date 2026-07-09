import type { Config } from 'tailwindcss';

export default {
  content: ['./index.html', './src/**/*.{ts,tsx}'],
  darkMode: 'class',
  theme: {
    extend: {
      colors: {
        shell: '#08111f',
        panel: '#101a2b',
        panel2: '#152238',
        line: '#263244',
        text: '#e5edf8',
        muted: '#91a1b8',
        aqua: '#5eead4',
        sky: '#38bdf8',
        magenta: '#f472b6',
        amber: '#facc15',
        coral: '#fb7185',
        leaf: '#86efac'
      },
      boxShadow: {
        tool: '0 18px 54px rgba(0, 0, 0, 0.28)',
        focus: '0 0 0 1px rgba(94, 234, 212, 0.28), 0 20px 64px rgba(6, 15, 28, 0.36)'
      }
    }
  },
  plugins: []
} satisfies Config;
