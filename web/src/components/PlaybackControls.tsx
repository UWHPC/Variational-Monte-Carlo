import type { CSSProperties } from 'react';
import type { PlaybackSpeed } from '../types/simulation';

interface PlaybackControlsProps {
  totalFrames: number;
  currentFrameIndex: number;
  isPlaying: boolean;
  speed: PlaybackSpeed;
  className?: string;
  compact?: boolean;
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
  className,
  compact = false,
  onPlay,
  onPause,
  onPrevious,
  onNext,
  onSeek,
  onSpeedChange,
}: PlaybackControlsProps) {
  const canStep = totalFrames > 1;
  const maxFrameIndex = Math.max(totalFrames - 1, 0);
  const sliderPercent =
    maxFrameIndex > 0 ? (currentFrameIndex / maxFrameIndex) * 100 : 0;
  const sliderStyle = {
    '--slider-progress': `${sliderPercent}%`,
  } as CSSProperties;

  return (
    <section
      className={`panel-card ${compact ? 'playback-compact' : ''} ${className ?? ''}`.trim()}
    >
      <h2 className="panel-title">Playback</h2>

      <div className="controls-row">
        <button
          className="control-btn"
          type="button"
          onClick={onPrevious}
          disabled={!canStep}
        >
          {compact ? 'Prev' : 'Previous'}
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
          className="playback-slider"
          type="range"
          min={0}
          max={maxFrameIndex}
          step={1}
          value={currentFrameIndex}
          style={sliderStyle}
          onChange={(event) => onSeek(Number(event.target.value))}
          disabled={!canStep}
        />
        <div className="range-labels">
          <span>0</span>
          <span>{maxFrameIndex}</span>
        </div>
      </div>

      <div className="speed-row">
        {!compact ? <span className="label-muted">Speed</span> : null}
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
