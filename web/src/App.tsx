import { CircleDot, Dna, Moon, Sun } from 'lucide-react';
import { useEffect, useState } from 'react';
import { Link, Route, Routes } from 'react-router-dom';
import Home from './routes/Home';
import ResultPage from './routes/ResultPage';
import UploadPage from './routes/UploadPage';

export default function App() {
  const [light, setLight] = useState(false);

  useEffect(() => {
    document.documentElement.classList.toggle('light', light);
  }, [light]);

  return (
    <div className="min-h-screen bg-shell text-text">
      <header className="sticky top-0 z-20 border-b border-line bg-shell/88 backdrop-blur-xl">
        <div className="mx-auto flex max-w-[1540px] items-center justify-between gap-4 px-5 py-3">
          <Link to="/" className="flex items-center gap-3">
            <span className="grid h-9 w-9 place-items-center rounded-md border border-aqua/40 bg-panel2 shadow-focus">
              <Dna className="h-5 w-5 text-aqua" aria-hidden />
            </span>
            <div>
              <div className="text-base font-semibold tracking-normal">Mito-Architect</div>
              <div className="flex items-center gap-1.5 text-xs text-muted">
                <CircleDot className="h-3 w-3 text-leaf" aria-hidden />
                mtDNA analysis console
              </div>
            </div>
          </Link>
          <button
            type="button"
            onClick={() => setLight((value) => !value)}
            className="inline-flex h-9 w-9 items-center justify-center rounded-md border border-line bg-panel2 text-muted hover:border-aqua hover:text-text"
            title="Toggle theme"
          >
            {light ? <Moon className="h-4 w-4" /> : <Sun className="h-4 w-4" />}
          </button>
        </div>
      </header>
      <Routes>
        <Route path="/" element={<Home />} />
        <Route path="/upload" element={<UploadPage />} />
        <Route path="/result/:jobId" element={<ResultPage />} />
      </Routes>
    </div>
  );
}
