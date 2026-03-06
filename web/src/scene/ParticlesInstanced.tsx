import { useLayoutEffect, useMemo, useRef } from 'react';
import * as THREE from 'three';
import type { InstancedMesh } from 'three';

interface ParticlesInstancedProps {
  numParticles: number;
  boxLength: number;
  positions: ArrayLike<number>;
}

export function ParticlesInstanced({
  numParticles,
  boxLength,
  positions,
}: ParticlesInstancedProps) {
  const meshRef = useRef<InstancedMesh>(null);
  const radius = useMemo(() => Math.max(0.08, boxLength * 0.03), [boxLength]);

  useLayoutEffect(() => {
    const mesh = meshRef.current;
    if (!mesh) {
      return;
    }
    mesh.instanceMatrix.setUsage(THREE.DynamicDrawUsage);
  }, []);

  useLayoutEffect(() => {
    const mesh = meshRef.current;
    if (!mesh) {
      return;
    }

    const expectedLength = numParticles * 3;
    if (positions.length !== expectedLength) {
      return;
    }

    const halfLength = boxLength / 2;
    const matrixArray = mesh.instanceMatrix.array as Float32Array;

    for (let i = 0; i < numParticles; i += 1) {
      const index = i * 3;
      const matrixOffset = i * 16;
      matrixArray[matrixOffset] = 1;
      matrixArray[matrixOffset + 1] = 0;
      matrixArray[matrixOffset + 2] = 0;
      matrixArray[matrixOffset + 3] = 0;
      matrixArray[matrixOffset + 4] = 0;
      matrixArray[matrixOffset + 5] = 1;
      matrixArray[matrixOffset + 6] = 0;
      matrixArray[matrixOffset + 7] = 0;
      matrixArray[matrixOffset + 8] = 0;
      matrixArray[matrixOffset + 9] = 0;
      matrixArray[matrixOffset + 10] = 1;
      matrixArray[matrixOffset + 11] = 0;
      matrixArray[matrixOffset + 12] = positions[index] - halfLength;
      matrixArray[matrixOffset + 13] = positions[index + 1] - halfLength;
      matrixArray[matrixOffset + 14] = positions[index + 2] - halfLength;
      matrixArray[matrixOffset + 15] = 1;
    }

    mesh.count = numParticles;
    mesh.instanceMatrix.needsUpdate = true;
  }, [boxLength, numParticles, positions]);

  return (
    <instancedMesh
      ref={meshRef}
      args={[undefined, undefined, numParticles]}
      frustumCulled={false}
    >
      <sphereGeometry args={[radius, 18, 18]} />
      <meshStandardMaterial color="#57b8ff" roughness={0.35} metalness={0.08} />
    </instancedMesh>
  );
}
