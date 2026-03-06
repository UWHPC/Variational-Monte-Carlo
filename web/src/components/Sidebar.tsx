import type {
  ReplayData,
  SimulationFrame,
} from '../types/simulation';
import type { PlaybackState } from '../hooks/usePlayback';
import { EnergyChart } from './EnergyChart';
import { PlaybackControls } from './PlaybackControls';
import { StatsPanel } from './StatsPanel';

interface SidebarProps {
  replay: ReplayData;
  currentFrame: SimulationFrame;
  playback: PlaybackState;
}

export function Sidebar({ replay, currentFrame, playback }: SidebarProps) {
  return (
    <aside className="sidebar">
      <StatsPanel
        init={replay.init}
        currentFrame={currentFrame}
        currentFrameIndex={playback.currentFrameIndex}
        totalFrames={replay.frames.length}
        done={replay.done}
      />
      <PlaybackControls
        totalFrames={replay.frames.length}
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
