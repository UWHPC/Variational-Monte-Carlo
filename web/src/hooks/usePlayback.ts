import { useCallback, useEffect, useState } from 'react';
import type { PlaybackSpeed } from '../types/simulation';

export interface PlaybackState {
  currentFrameIndex: number;
  isPlaying: boolean;
  speed: PlaybackSpeed;
  setSpeed: (speed: PlaybackSpeed) => void;
  play: () => void;
  pause: () => void;
  previous: () => void;
  next: () => void;
  seek: (index: number) => void;
  lastFrameIndex: number;
}

const BASE_FRAMES_PER_SECOND = 12;

function clamp(value: number, min: number, max: number): number {
  return Math.min(max, Math.max(min, value));
}

export function usePlayback(frameCount: number): PlaybackState {
  const lastFrameIndex = Math.max(0, frameCount - 1);

  const [currentFrameIndex, setCurrentFrameIndex] = useState(0);
  const [isPlaying, setIsPlaying] = useState(false);
  const [speed, setSpeed] = useState<PlaybackSpeed>(1);

  useEffect(() => {
    if (frameCount === 0) {
      setCurrentFrameIndex(0);
      setIsPlaying(false);
      return;
    }

    setCurrentFrameIndex((index) => clamp(index, 0, lastFrameIndex));
  }, [frameCount, lastFrameIndex]);

  useEffect(() => {
    if (isPlaying && currentFrameIndex >= lastFrameIndex) {
      setIsPlaying(false);
    }
  }, [currentFrameIndex, isPlaying, lastFrameIndex]);

  useEffect(() => {
    if (!isPlaying || frameCount <= 1) {
      return;
    }

    const intervalMs = Math.max(20, 1000 / (BASE_FRAMES_PER_SECOND * speed));
    const intervalId = window.setInterval(() => {
      setCurrentFrameIndex((index) => Math.min(index + 1, lastFrameIndex));
    }, intervalMs);

    return () => {
      window.clearInterval(intervalId);
    };
  }, [frameCount, isPlaying, lastFrameIndex, speed]);

  const play = useCallback(() => {
    if (frameCount <= 1) {
      setIsPlaying(false);
      return;
    }

    setCurrentFrameIndex((index) => (index >= lastFrameIndex ? 0 : index));
    setIsPlaying(true);
  }, [frameCount, lastFrameIndex]);

  const pause = useCallback(() => {
    setIsPlaying(false);
  }, []);

  const previous = useCallback(() => {
    setIsPlaying(false);
    setCurrentFrameIndex((index) => Math.max(index - 1, 0));
  }, []);

  const next = useCallback(() => {
    setIsPlaying(false);
    setCurrentFrameIndex((index) => Math.min(index + 1, lastFrameIndex));
  }, [lastFrameIndex]);

  const seek = useCallback(
    (index: number) => {
      setIsPlaying(false);
      setCurrentFrameIndex(clamp(Math.round(index), 0, lastFrameIndex));
    },
    [lastFrameIndex],
  );

  return {
    currentFrameIndex,
    isPlaying,
    speed,
    setSpeed,
    play,
    pause,
    previous,
    next,
    seek,
    lastFrameIndex,
  };
}
