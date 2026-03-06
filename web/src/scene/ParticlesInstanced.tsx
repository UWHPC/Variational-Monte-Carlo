import { useLayoutEffect, useMemo, useRef } from 'react';
import * as THREE from 'three';
import type { InstancedMesh } from 'three';

interface ParticlesInstancedProps {
  numParticles: number;
  boxLength: number;
  positions: number[];
}

export function ParticlesInstanced({
  numParticles,
  boxLength,
  positions,
}: ParticlesInstancedProps) {
  const meshRef = useRef<InstancedMesh>(null);
  const tempObject = useMemo(() => new THREE.Object3D(), []);
  const radius = useMemo(() => Math.max(0.08, boxLength * 0.03), [boxLength]);

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

    for (let i = 0; i < numParticles; i += 1) {
      const index = i * 3;
      tempObject.position.set(
        positions[index] - halfLength,
        positions[index + 1] - halfLength,
        positions[index + 2] - halfLength,
      );
      tempObject.updateMatrix();
      mesh.setMatrixAt(i, tempObject.matrix);
    }

    mesh.instanceMatrix.needsUpdate = true;
  }, [boxLength, numParticles, positions, tempObject]);

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
