import type {
  SimulationDone,
  SimulationFrame,
  SimulationInit,
} from '../types/simulation';
import {
  formatInteger,
  formatNumber,
  formatOptionalNumber,
  formatPercent,
} from '../lib/format';

interface StatsPanelProps {
  init: SimulationInit;
  currentFrame: SimulationFrame;
  currentFrameIndex: number;
  totalFrames: number;
  done?: SimulationDone;
}

interface MetricRowProps {
  label: string;
  value: string;
}

function MetricRow({ label, value }: MetricRowProps) {
  return (
    <div className="metric-row">
      <span className="metric-label">{label}</span>
      <span className="metric-value">{value}</span>
    </div>
  );
}

export function StatsPanel({
  init,
  currentFrame,
  currentFrameIndex,
  totalFrames,
  done,
}: StatsPanelProps) {
  return (
    <section className="panel-card">
      <h2 className="panel-title">Simulation Stats</h2>

      <div className="metric-grid">
        <MetricRow
          label="Frame"
          value={`${formatInteger(currentFrameIndex)} / ${formatInteger(
            Math.max(totalFrames - 1, 0),
          )}`}
        />
        <MetricRow label="Step" value={formatInteger(currentFrame.step)} />
        <MetricRow
          label="Local Energy"
          value={formatNumber(currentFrame.localEnergy, 6)}
        />
        <MetricRow
          label="Running Mean"
          value={formatNumber(currentFrame.meanEnergy, 6)}
        />
        <MetricRow
          label="Std Error"
          value={formatOptionalNumber(
            currentFrame.standardErrorAvailable,
            currentFrame.standardError,
            6,
          )}
        />
        <MetricRow
          label="Acceptance"
          value={formatPercent(currentFrame.acceptanceRate, 2)}
        />
      </div>

      <div className="panel-divider" />

      <div className="metric-grid">
        <MetricRow label="Run ID" value={init.runId} />
        <MetricRow
          label="Particles"
          value={formatInteger(init.numParticles)}
        />
        <MetricRow
          label="Box Length"
          value={formatNumber(init.boxLength, 3)}
        />
        <MetricRow
          label="Warmup Steps"
          value={formatInteger(init.warmupSteps)}
        />
        <MetricRow
          label="Measure Steps"
          value={formatInteger(init.measureSteps)}
        />
        <MetricRow
          label="Step Size"
          value={formatNumber(init.stepSize, 4)}
        />
        <MetricRow label="Block Size" value={formatInteger(init.blockSize)} />
        <MetricRow label="Seed" value={formatInteger(init.seed)} />
      </div>

      {done ? (
        <>
          <div className="panel-divider" />
          <div className="metric-grid">
            <MetricRow
              label="Final Mean"
              value={formatNumber(done.finalMeanEnergy, 6)}
            />
            <MetricRow
              label="Final Std Error"
              value={formatOptionalNumber(
                done.finalStandardErrorAvailable,
                done.finalStandardError,
                6,
              )}
            />
            <MetricRow
              label="Final Acceptance"
              value={formatPercent(done.finalAcceptanceRate, 2)}
            />
          </div>
        </>
      ) : null}
    </section>
  );
}
