import type {
  LoaderResult,
  ReplayData,
  SimulationDone,
  SimulationFrame,
  SimulationInit,
} from '../types/simulation';

type JsonRecord = Record<string, unknown>;

function isRecord(value: unknown): value is JsonRecord {
  return typeof value === 'object' && value !== null;
}

function isFiniteNumber(value: unknown): value is number {
  return typeof value === 'number' && Number.isFinite(value);
}

function readString(record: JsonRecord, key: string, line: number): string {
  const value = record[key];
  if (typeof value !== 'string' || value.length === 0) {
    throw new Error(`Line ${line}: "${key}" must be a non-empty string.`);
  }
  return value;
}

function readBoolean(record: JsonRecord, key: string, line: number): boolean {
  const value = record[key];
  if (typeof value !== 'boolean') {
    throw new Error(`Line ${line}: "${key}" must be a boolean.`);
  }
  return value;
}

function readNumber(record: JsonRecord, key: string, line: number): number {
  const value = record[key];
  if (!isFiniteNumber(value)) {
    throw new Error(`Line ${line}: "${key}" must be a finite number.`);
  }
  return value;
}

function readNumberArray(record: JsonRecord, key: string, line: number): number[] {
  const value = record[key];
  if (!Array.isArray(value)) {
    throw new Error(`Line ${line}: "${key}" must be an array.`);
  }

  const parsed: number[] = [];
  for (let i = 0; i < value.length; i += 1) {
    const element = value[i];
    if (!isFiniteNumber(element)) {
      throw new Error(
        `Line ${line}: "${key}[${i}]" must be a finite number.`,
      );
    }
    parsed.push(element);
  }

  return parsed;
}

function parseInit(record: JsonRecord, line: number): SimulationInit {
  const init: SimulationInit = {
    type: 'init',
    runId: readString(record, 'runId', line),
    numParticles: readNumber(record, 'numParticles', line),
    boxLength: readNumber(record, 'boxLength', line),
    warmupSteps: readNumber(record, 'warmupSteps', line),
    measureSteps: readNumber(record, 'measureSteps', line),
    stepSize: readNumber(record, 'stepSize', line),
    blockSize: readNumber(record, 'blockSize', line),
    seed: readNumber(record, 'seed', line),
  };

  if (init.numParticles <= 0) {
    throw new Error(`Line ${line}: "numParticles" must be > 0.`);
  }
  if (init.boxLength <= 0) {
    throw new Error(`Line ${line}: "boxLength" must be > 0.`);
  }

  return init;
}

function parseFrame(
  record: JsonRecord,
  line: number,
  init: SimulationInit,
): SimulationFrame {
  const standardErrorAvailable = readBoolean(
    record,
    'standardErrorAvailable',
    line,
  );
  let standardError = 0;

  if (standardErrorAvailable) {
    standardError = readNumber(record, 'standardError', line);
  } else if (
    Object.prototype.hasOwnProperty.call(record, 'standardError') &&
    !isFiniteNumber(record.standardError)
  ) {
    throw new Error(
      `Line ${line}: "standardError" must be a finite number when provided.`,
    );
  }

  const positions = readNumberArray(record, 'positions', line);
  const expected = init.numParticles * 3;

  if (positions.length !== expected) {
    throw new Error(
      `Line ${line}: "positions" length must be ${expected}, got ${positions.length}.`,
    );
  }

  return {
    type: 'frame',
    step: readNumber(record, 'step', line),
    accepted: readNumber(record, 'accepted', line),
    proposed: readNumber(record, 'proposed', line),
    acceptanceRate: readNumber(record, 'acceptanceRate', line),
    localEnergy: readNumber(record, 'localEnergy', line),
    meanEnergy: readNumber(record, 'meanEnergy', line),
    standardErrorAvailable,
    standardError,
    positions,
  };
}

function parseDone(record: JsonRecord, line: number): SimulationDone {
  const finalStandardErrorAvailable = readBoolean(
    record,
    'finalStandardErrorAvailable',
    line,
  );

  let finalStandardError = 0;
  if (finalStandardErrorAvailable) {
    finalStandardError = readNumber(record, 'finalStandardError', line);
  } else if (
    Object.prototype.hasOwnProperty.call(record, 'finalStandardError') &&
    !isFiniteNumber(record.finalStandardError)
  ) {
    throw new Error(
      `Line ${line}: "finalStandardError" must be a finite number when provided.`,
    );
  }

  return {
    type: 'done',
    totalAccepted: readNumber(record, 'totalAccepted', line),
    totalProposed: readNumber(record, 'totalProposed', line),
    finalAcceptanceRate: readNumber(record, 'finalAcceptanceRate', line),
    finalMeanEnergy: readNumber(record, 'finalMeanEnergy', line),
    finalStandardErrorAvailable,
    finalStandardError,
  };
}

export function parseReplayJsonl(rawText: string): LoaderResult {
  const lines = rawText.split(/\r?\n/);

  let init: SimulationInit | undefined;
  const frames: SimulationFrame[] = [];
  let done: SimulationDone | undefined;

  for (let lineIndex = 0; lineIndex < lines.length; lineIndex += 1) {
    const line = lines[lineIndex].trim();
    const lineNumber = lineIndex + 1;

    if (!line) {
      continue;
    }

    let parsed: unknown;
    try {
      parsed = JSON.parse(line);
    } catch {
      return {
        ok: false,
        error: `Line ${lineNumber}: invalid JSON.`,
      };
    }

    if (!isRecord(parsed)) {
      return {
        ok: false,
        error: `Line ${lineNumber}: each message must be a JSON object.`,
      };
    }

    const type = parsed.type;
    if (type === 'init') {
      if (init) {
        return {
          ok: false,
          error: `Line ${lineNumber}: multiple init messages found; exactly one is required.`,
        };
      }

      try {
        init = parseInit(parsed, lineNumber);
      } catch (error) {
        return {
          ok: false,
          error:
            error instanceof Error
              ? error.message
              : `Line ${lineNumber}: invalid init message.`,
        };
      }
      continue;
    }

    if (type === 'frame') {
      if (!init) {
        return {
          ok: false,
          error: `Line ${lineNumber}: frame message encountered before init.`,
        };
      }

      try {
        frames.push(parseFrame(parsed, lineNumber, init));
      } catch (error) {
        return {
          ok: false,
          error:
            error instanceof Error
              ? error.message
              : `Line ${lineNumber}: invalid frame message.`,
        };
      }
      continue;
    }

    if (type === 'done') {
      if (!init) {
        return {
          ok: false,
          error: `Line ${lineNumber}: done message encountered before init.`,
        };
      }
      if (done) {
        return {
          ok: false,
          error: `Line ${lineNumber}: multiple done messages found; at most one is allowed.`,
        };
      }

      try {
        done = parseDone(parsed, lineNumber);
      } catch (error) {
        return {
          ok: false,
          error:
            error instanceof Error
              ? error.message
              : `Line ${lineNumber}: invalid done message.`,
        };
      }
      continue;
    }

    return {
      ok: false,
      error: `Line ${lineNumber}: unknown message type "${String(type)}".`,
    };
  }

  if (!init) {
    return {
      ok: false,
      error: 'Replay is invalid: exactly one init message is required.',
    };
  }

  if (frames.length === 0) {
    return {
      ok: false,
      error: 'Replay is invalid: at least one frame message is required.',
    };
  }

  const replayData: ReplayData = {
    init,
    frames,
    done,
  };

  return {
    ok: true,
    data: replayData,
  };
}

export async function loadFrames(url: string): Promise<LoaderResult> {
  let response: Response;

  try {
    response = await fetch(url);
  } catch (error) {
    return {
      ok: false,
      error: `Failed to fetch replay file from "${url}": ${
        error instanceof Error ? error.message : 'network error'
      }.`,
    };
  }

  if (!response.ok) {
    return {
      ok: false,
      error: `Failed to fetch replay file from "${url}": HTTP ${response.status}.`,
    };
  }

  let rawText = '';
  try {
    rawText = await response.text();
  } catch (error) {
    return {
      ok: false,
      error: `Failed reading replay response: ${
        error instanceof Error ? error.message : 'unknown error'
      }.`,
    };
  }

  return parseReplayJsonl(rawText);
}
