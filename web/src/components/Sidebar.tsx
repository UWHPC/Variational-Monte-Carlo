import type {
  ReplayData,
  SimulationFrame,
} from '../types/simulation';
import type { PlaybackState } from '../hooks/usePlayback';
import { formatInteger } from '../lib/format';
import { EnergyChart } from './EnergyChart';
import { PlaybackControls } from './PlaybackControls';
import { StatsPanel } from './StatsPanel';

interface SidebarProps {
  replay: ReplayData;
  currentFrame: SimulationFrame;
  playback: PlaybackState;
}

export function Sidebar({ replay, currentFrame, playback }: SidebarProps) {
  const playState = playback.isPlaying ? 'Playing' : 'Paused';
  const totalFrames = replay.frames.length;
  const frameLabel = `${formatInteger(playback.currentFrameIndex + 1)} / ${formatInteger(totalFrames)}`;

  return (
    <aside className="sidebar">
      <section className="panel-card status-strip">
        <div className="status-chip-group">
          <div className="status-chip">
            <span className="status-label">Run</span>
            <span className="status-value">{replay.init.runId}</span>
          </div>
          <div className="status-chip">
            <span className="status-label">Frame</span>
            <span className="status-value mono">{frameLabel}</span>
          </div>
          <div className="status-chip">
            <span className="status-label">Speed</span>
            <span className="status-value mono">{playback.speed}x</span>
          </div>
          <div
            className={`status-chip ${playback.isPlaying ? 'status-chip-live' : ''}`}
          >
            <span className="status-label">State</span>
            <span className="status-value">{playState}</span>
          </div>
        </div>
      </section>

      <StatsPanel
        init={replay.init}
        currentFrame={currentFrame}
        currentFrameIndex={playback.currentFrameIndex}
        totalFrames={totalFrames}
        done={replay.done}
      />
      <PlaybackControls
        totalFrames={totalFrames}
        currentFrameIndex={playback.currentFrameIndex}
        isPlaying={playback.isPlaying}
        speed={playback.speed}
        onPlay={playback.play}
        onPause={playback.pause}
        onPrevious={playback.previous}
        onNext={playback.next}
        onSeek={playback.seek}
        onSpeedChange={playback.setSpeed}
      />
      <EnergyChart
        frames={replay.frames}
        currentFrameIndex={playback.currentFrameIndex}
      />
    </aside>
  );
}
