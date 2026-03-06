import { useMemo } from 'react';
import { Canvas } from '@react-three/fiber';
import { OrbitControls } from '@react-three/drei';
import { BoxWireframe } from './BoxWireframe';
import { ParticlesInstanced } from './ParticlesInstanced';

interface SimulationCanvasProps {
  boxLength: number;
  numParticles: number;
  positions: ArrayLike<number>;
}

function getDefaultCameraPosition(boxLength: number): [number, number, number] {
  const distance = Math.max(8, boxLength * 1.8);
  return [distance, distance * 0.75, distance * 0.95];
}

export function SimulationCanvas({
  boxLength,
  numParticles,
  positions,
}: SimulationCanvasProps) {
  const cameraPosition = useMemo(
    () => getDefaultCameraPosition(boxLength),
    [boxLength],
  );

  return (
    <div className="canvas-shell">
      <Canvas
        dpr={[1, 2]}
        camera={{
          position: cameraPosition,
          fov: 48,
          near: 0.1,
          far: Math.max(200, boxLength * 20),
        }}
        gl={{ antialias: true }}
      >
        <color attach="background" args={['#0b0d11']} />
        <ambientLight intensity={0.46} />
        <directionalLight
          position={[boxLength * 1.2, boxLength * 1.8, boxLength * 0.9]}
          intensity={1.1}
        />
        <BoxWireframe boxLength={boxLength} />
        <ParticlesInstanced
          numParticles={numParticles}
          boxLength={boxLength}
          positions={positions}
        />
        <OrbitControls
          makeDefault
          enableDamping
          dampingFactor={0.08}
          minDistance={boxLength * 0.4}
          maxDistance={boxLength * 6}
        />
      </Canvas>
    </div>
  );
}
