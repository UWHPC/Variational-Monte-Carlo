import type { PlaybackSpeed } from '../types/simulation';

interface PlaybackControlsProps {
  totalFrames: number;
  currentFrameIndex: number;
  isPlaying: boolean;
  speed: PlaybackSpeed;
  onPlay: () => void;
  onPause: () => void;
  onPrevious: () => void;
  onNext: () => void;
  onSeek: (index: number) => void;
  onSpeedChange: (speed: PlaybackSpeed) => void;
}

const SPEED_OPTIONS: PlaybackSpeed[] = [0.25, 1, 4, 10];

export function PlaybackControls({
  totalFrames,
  currentFrameIndex,
  isPlaying,
  speed,
  onPlay,
  onPause,
  onPrevious,
  onNext,
  onSeek,
  onSpeedChange,
}: PlaybackControlsProps) {
  const canStep = totalFrames > 1;

  return (
    <section className="panel-card">
      <h2 className="panel-title">Playback</h2>

      <div className="controls-row">
        <button
          className="control-btn"
          type="button"
          onClick={onPrevious}
          disabled={!canStep}
        >
          Previous
        </button>
        {isPlaying ? (
          <button
            className="control-btn"
            type="button"
            onClick={onPause}
            disabled={!canStep}
          >
            Pause
          </button>
        ) : (
          <button
            className="control-btn"
            type="button"
            onClick={onPlay}
            disabled={!canStep}
          >
            Play
          </button>
        )}
        <button
          className="control-btn"
          type="button"
          onClick={onNext}
          disabled={!canStep}
        >
          Next
        </button>
      </div>

      <div className="slider-wrap">
        <input
          type="range"
          min={0}
          max={Math.max(totalFrames - 1, 0)}
          step={1}
          value={currentFrameIndex}
          onChange={(event) => onSeek(Number(event.target.value))}
          disabled={!canStep}
        />
        <div className="range-labels">
          <span>0</span>
          <span>{Math.max(totalFrames - 1, 0)}</span>
        </div>
      </div>

      <div className="speed-row">
        <span className="label-muted">Speed</span>
        <select
          value={speed}
          onChange={(event) =>
            onSpeedChange(Number(event.target.value) as PlaybackSpeed)
          }
        >
          {SPEED_OPTIONS.map((option) => (
            <option key={option} value={option}>
              {option}x
            </option>
          ))}
        </select>
      </div>
    </section>
  );
}
