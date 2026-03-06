export interface SimulationInit {
  type: 'init';
  runId: string;
  numParticles: number;
  boxLength: number;
  warmupSteps: number;
  measureSteps: number;
  stepSize: number;
  blockSize: number;
  seed: number;
}

export interface SimulationFrame {
  type: 'frame';
  step: number;
  accepted: number;
  proposed: number;
  acceptanceRate: number;
  localEnergy: number;
  meanEnergy: number;
  standardErrorAvailable: boolean;
  standardError: number;
  positions: number[];
}

export interface SimulationDone {
  type: 'done';
  totalAccepted: number;
  totalProposed: number;
  finalAcceptanceRate: number;
  finalMeanEnergy: number;
  finalStandardErrorAvailable: boolean;
  finalStandardError: number;
}

export interface ReplayData {
  init: SimulationInit;
  frames: SimulationFrame[];
  done?: SimulationDone;
}

export type PlaybackSpeed = 0.25 | 1 | 4 | 10;

export type LoaderResult =
  | {
      ok: true;
      data: ReplayData;
    }
  | {
      ok: false;
      error: string;
    };
