import { useCallback, useEffect, useMemo, useState } from 'react';
import { ErrorState } from './components/ErrorState';
import { LoadingState } from './components/LoadingState';
import { PlaybackControls } from './components/PlaybackControls';
import { Sidebar } from './components/Sidebar';
import { usePlayback } from './hooks/usePlayback';
import { loadFrames } from './lib/loadFrames';
import { SimulationCanvas } from './scene/SimulationCanvas';
import type { ReplayData } from './types/simulation';

type LoadStatus =
  | { state: 'loading' }
  | { state: 'error'; message: string }
  | { state: 'ready'; replay: ReplayData };

export default function App() {
  const [loadStatus, setLoadStatus] = useState<LoadStatus>({ state: 'loading' });

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

    const result = await loadFrames('/sample-run.jsonl');
    if (!result.ok) {
      setLoadStatus({ state: 'error', message: result.error });
      return;
    }

    setLoadStatus({ state: 'ready', replay: result.data });
  }, []);

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
    return <LoadingState />;
  }

  if (loadStatus.state === 'error') {
    return <ErrorState message={loadStatus.message} onRetry={loadReplay} />;
  }

  if (!currentFrame) {
    return (
      <ErrorState
        message="Replay loaded but no current frame is available."
        onRetry={loadReplay}
      />
    );
  }

  const replay = loadStatus.replay;

  return (
    <div className="app-root">
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
