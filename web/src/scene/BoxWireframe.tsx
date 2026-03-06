interface BoxWireframeProps {
  boxLength: number;
}

export function BoxWireframe({ boxLength }: BoxWireframeProps) {
  return (
    <mesh>
      <boxGeometry args={[boxLength, boxLength, boxLength]} />
      <meshBasicMaterial
        color="#2b4059"
        wireframe
        transparent
        opacity={0.9}
      />
    </mesh>
  );
}
