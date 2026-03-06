const LOCALE = 'en-US';

export function formatInteger(value: number): string {
  if (!Number.isFinite(value)) {
    return 'NaN';
  }
  return Math.trunc(value).toLocaleString(LOCALE);
}

export function formatNumber(value: number, fractionDigits = 4): string {
  if (!Number.isFinite(value)) {
    return 'NaN';
  }
  return value.toLocaleString(LOCALE, {
    minimumFractionDigits: fractionDigits,
    maximumFractionDigits: fractionDigits,
  });
}

export function formatPercent(value: number, fractionDigits = 2): string {
  if (!Number.isFinite(value)) {
    return 'NaN';
  }
  return `${(value * 100).toFixed(fractionDigits)}%`;
}

export function formatOptionalNumber(
  available: boolean,
  value: number,
  fractionDigits = 4,
  unavailableText = 'unavailable',
): string {
  if (!available) {
    return unavailableText;
  }
  return formatNumber(value, fractionDigits);
}
