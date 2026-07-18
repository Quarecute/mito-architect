export type AnalysisTabNavigationKey = 'ArrowLeft' | 'ArrowRight' | 'Home' | 'End';

export function nextTabIndex(
  key: string,
  currentIndex: number,
  tabCount: number
): number | undefined {
  if (!Number.isInteger(currentIndex) || !Number.isInteger(tabCount) || tabCount <= 0) {
    return undefined;
  }
  if (key === 'Home') return 0;
  if (key === 'End') return tabCount - 1;
  if (key === 'ArrowRight') return (currentIndex + 1) % tabCount;
  if (key === 'ArrowLeft') return (currentIndex - 1 + tabCount) % tabCount;
  return undefined;
}

export function normalizeExactEvidenceFilter(value: string | undefined): string | undefined {
  const normalized = value?.trim();
  return normalized ? normalized : undefined;
}
