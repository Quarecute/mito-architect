import assert from 'node:assert/strict';
import test from 'node:test';

import {
  normalizePdbId,
  rcsbCandidates
} from './proteinStructureSources.ts';

test('RCSB sources are deterministic and local-first', () => {
  const candidates = rcsbCandidates('8h9s', 'N');
  assert.equal(candidates.length, 3);
  assert.equal(candidates[0].url, '/structures/8h9s.bcif');
  assert.equal(candidates[0].chain, 'N');
  assert.equal(candidates[0].confidenceSemantics, 'bfactor');
  assert.equal(candidates[1].url, 'https://models.rcsb.org/8h9s.bcif');
  assert.equal(candidates[2].url, 'https://files.rcsb.org/download/8H9S.cif');
});

test('structure identifiers are validated and normalized', () => {
  assert.equal(normalizePdbId(' 9i4i '), '9I4I');
  assert.equal(normalizePdbId('Q00360'), undefined);
  assert.equal(normalizePdbId('AF-P00846-F1'), undefined);
});
