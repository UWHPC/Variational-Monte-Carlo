import { useCallback, useEffect, useMemo, useState, type ChangeEvent } from 'react';
import { ErrorState } from './components/ErrorState';
import { LoadingState } from './components/LoadingState';
import { PlaybackControls } from './components/PlaybackControls';
import { Sidebar } from './components/Sidebar';
import { usePlayback } from './hooks/usePlayback';
import { loadFrames, parseReplayJsonl } from './lib/loadFrames';
import { SimulationCanvas } from './scene/SimulationCanvas';
import type { ReplayData } from './types/simulation';

type LoadStatus =
  | { state: 'loading' }
  | { state: 'error'; message: string }
  | { state: 'ready'; replay: ReplayData };

interface ReplaySourceControlsProps {
  sourceLabel: string;
  isBusy: boolean;
  onUpload: (event: ChangeEvent<HTMLInputElement>) => void;
  onLoadSample: () => void;
}

function ReplaySourceControls({
  sourceLabel,
  isBusy,
  onUpload,
  onLoadSample,
}: ReplaySourceControlsProps) {
  return (
    <section className="panel-card replay-toolbar">
      <div className="replay-toolbar-actions">
        <label className="toolbar-btn" role="button">
          Upload JSONL
          <input
            className="replay-upload-input"
            type="file"
            accept=".jsonl,application/x-ndjson,application/json,text/plain"
            onChange={onUpload}
            disabled={isBusy}
          />
        </label>
        <button
          className="toolbar-btn"
          type="button"
          onClick={onLoadSample}
          disabled={isBusy}
        >
          Load Sample
        </button>
      </div>
      <div className="replay-source-label">{sourceLabel}</div>
    </section>
  );
}

export default function App() {
  const samplePath = '/sample-run.jsonl';
  const [loadStatus, setLoadStatus] = useState<LoadStatus>({ state: 'loading' });
  const [sourceLabel, setSourceLabel] = useState<string>(`Sample: ${samplePath}`);

  const frameCount =
    loadStatus.state === 'ready' ? loadStatus.replay.frames.length : 0;
  const playback = usePlayback(frameCount);

  const currentFrame = useMemo(() => {
    if (loadStatus.state !== 'ready') {
      return undefined;
    }
    return (
      loadStatus.replay.frames[playback.currentFrameIndex] ??
      loadStatus.replay.frames[0]
    );
  }, [loadStatus, playback.currentFrameIndex]);

  const loadReplay = useCallback(async () => {
    setLoadStatus({ state: 'loading' });

    const result = await loadFrames(samplePath);
    if (!result.ok) {
      setLoadStatus({ state: 'error', message: result.error });
      return;
    }

    setSourceLabel(`Sample: ${samplePath}`);
    setLoadStatus({ state: 'ready', replay: result.data });
  }, [samplePath]);

  const loadReplayFromUpload = useCallback(
    async (event: ChangeEvent<HTMLInputElement>) => {
      const input = event.currentTarget;
      const file = input.files?.[0];
      input.value = '';

      if (!file) {
        return;
      }

      setLoadStatus({ state: 'loading' });

      let rawText = '';
      try {
        rawText = await file.text();
      } catch (error) {
        setLoadStatus({
          state: 'error',
          message: `Failed to read "${file.name}": ${
            error instanceof Error ? error.message : 'unknown error'
          }.`,
        });
        return;
      }

      const result = parseReplayJsonl(rawText);
      if (!result.ok) {
        setLoadStatus({ state: 'error', message: result.error });
        return;
      }

      setSourceLabel(`File: ${file.name}`);
      setLoadStatus({ state: 'ready', replay: result.data });
    },
    [],
  );

  useEffect(() => {
    void loadReplay();
  }, [loadReplay]);

  useEffect(() => {
    const isTypingTarget = (target: EventTarget | null): boolean => {
      if (!(target instanceof HTMLElement)) {
        return false;
      }

      if (target.isContentEditable) {
        return true;
      }

      const tagName = target.tagName;
      return (
        tagName === 'INPUT' ||
        tagName === 'TEXTAREA' ||
        tagName === 'SELECT' ||
        tagName === 'BUTTON'
      );
    };

    const onKeyDown = (event: KeyboardEvent) => {
      if (isTypingTarget(event.target)) {
        return;
      }

      if (event.code === 'Space') {
        event.preventDefault();
        if (event.repeat) {
          return;
        }

        if (playback.isPlaying) {
          playback.pause();
        } else {
          playback.play();
        }
        return;
      }

      if (event.key === 'ArrowLeft') {
        event.preventDefault();
        playback.previous();
        return;
      }

      if (event.key === 'ArrowRight') {
        event.preventDefault();
        playback.next();
      }
    };

    window.addEventListener('keydown', onKeyDown);
    return () => {
      window.removeEventListener('keydown', onKeyDown);
    };
  }, [
    playback.isPlaying,
    playback.next,
    playback.pause,
    playback.play,
    playback.previous,
  ]);

  if (loadStatus.state === 'loading') {
    return (
      <div className="app-root">
        <ReplaySourceControls
          sourceLabel={sourceLabel}
          isBusy
          onUpload={loadReplayFromUpload}
          onLoadSample={loadReplay}
        />
        <div className="state-host">
          <LoadingState />
        </div>
      </div>
    );
  }

  if (loadStatus.state === 'error') {
    return (
      <div className="app-root">
        <ReplaySourceControls
          sourceLabel={sourceLabel}
          isBusy={false}
          onUpload={loadReplayFromUpload}
          onLoadSample={loadReplay}
        />
        <div className="state-host">
          <ErrorState message={loadStatus.message} onRetry={loadReplay} />
        </div>
      </div>
    );
  }

  if (!currentFrame) {
    return (
      <div className="app-root">
        <ReplaySourceControls
          sourceLabel={sourceLabel}
          isBusy={false}
          onUpload={loadReplayFromUpload}
          onLoadSample={loadReplay}
        />
        <div className="state-host">
          <ErrorState
            message="Replay loaded but no current frame is available."
            onRetry={loadReplay}
          />
        </div>
      </div>
    );
  }

  const replay = loadStatus.replay;

  return (
    <div className="app-root">
      <ReplaySourceControls
        sourceLabel={sourceLabel}
        isBusy={false}
        onUpload={loadReplayFromUpload}
        onLoadSample={loadReplay}
      />
      <main className="app-layout">
        <section className="scene-panel">
          <SimulationCanvas
            boxLength={replay.init.boxLength}
            numParticles={replay.init.numParticles}
            positions={currentFrame.positions}
          />
          <PlaybackControls
            className="scene-playback-floating"
            compact
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
        </section>
        <Sidebar replay={replay} currentFrame={currentFrame} playback={playback} />
      </main>
    </div>
  );
}
