import assert from 'node:assert/strict';
import test from 'node:test';

import {
  nextTabIndex,
  normalizeExactEvidenceFilter
} from './analysisNavigation.ts';

test('analysis tabs implement wrapping arrow and absolute Home/End navigation', () => {
  assert.equal(nextTabIndex('ArrowRight', 5, 6), 0);
  assert.equal(nextTabIndex('ArrowLeft', 0, 6), 5);
  assert.equal(nextTabIndex('Home', 4, 6), 0);
  assert.equal(nextTabIndex('End', 1, 6), 5);
  assert.equal(nextTabIndex('Enter', 1, 6), undefined);
  assert.equal(nextTabIndex('ArrowRight', 0, 0), undefined);
});

test('global evidence filters preserve exact identifiers and reject whitespace-only input', () => {
  assert.equal(normalizeExactEvidenceFilter(' molecule-01 '), 'molecule-01');
  assert.equal(normalizeExactEvidenceFilter('sv:deletion:10-20'), 'sv:deletion:10-20');
  assert.equal(normalizeExactEvidenceFilter('   '), undefined);
  assert.equal(normalizeExactEvidenceFilter(undefined), undefined);
});
