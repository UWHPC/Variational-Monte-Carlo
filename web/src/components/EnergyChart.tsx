import { useMemo } from 'react';
import {
  CartesianGrid,
  Legend,
  Line,
  LineChart,
  ReferenceLine,
  ResponsiveContainer,
  Tooltip,
  XAxis,
  YAxis,
} from 'recharts';
import type { SimulationFrame } from '../types/simulation';
import { formatNumber } from '../lib/format';

interface EnergyChartProps {
  frames: SimulationFrame[];
  currentFrameIndex: number;
}

interface ChartPoint {
  step: number;
  localEnergy: number;
  meanEnergy: number;
}

export function EnergyChart({ frames, currentFrameIndex }: EnergyChartProps) {
  const data = useMemo<ChartPoint[]>(
    () =>
      frames.map((frame) => ({
        step: frame.step,
        localEnergy: frame.localEnergy,
        meanEnergy: frame.meanEnergy,
      })),
    [frames],
  );

  const currentStep = frames[currentFrameIndex]?.step;

  return (
    <section className="panel-card chart-panel">
      <h2 className="panel-title">Energy History</h2>
      <div className="chart-shell">
        <ResponsiveContainer width="100%" height="100%">
          <LineChart data={data} margin={{ top: 8, right: 12, left: 8, bottom: 8 }}>
            <CartesianGrid strokeDasharray="3 3" stroke="#26344a" />
            <XAxis
              dataKey="step"
              tick={{ fill: '#93a4bd', fontSize: 11 }}
              axisLine={{ stroke: '#334861' }}
              tickLine={{ stroke: '#334861' }}
            />
            <YAxis
              tick={{ fill: '#93a4bd', fontSize: 11 }}
              axisLine={{ stroke: '#334861' }}
              tickLine={{ stroke: '#334861' }}
              domain={['auto', 'auto']}
            />
            <Tooltip
              contentStyle={{
                backgroundColor: '#111a26',
                border: '1px solid #2b3f59',
                borderRadius: 8,
              }}
              formatter={(value: unknown, name: string) => {
                const numericValue =
                  typeof value === 'number' ? value : Number(value);
                const label =
                  name === 'localEnergy' ? 'Local Energy' : 'Running Mean';
                return [formatNumber(numericValue, 6), label];
              }}
              labelFormatter={(value: number | string) => `Step ${value}`}
            />
            <Legend
              formatter={(value: string) =>
                value === 'localEnergy' ? 'Local Energy' : 'Running Mean'
              }
              wrapperStyle={{ color: '#d3deed', fontSize: 12 }}
            />
            {Number.isFinite(currentStep) ? (
              <ReferenceLine
                x={currentStep}
                stroke="#ffcc66"
                strokeDasharray="4 3"
              />
            ) : null}
            <Line
              type="monotone"
              dataKey="localEnergy"
              stroke="#56b5ff"
              strokeWidth={2}
              dot={false}
              isAnimationActive={false}
            />
            <Line
              type="monotone"
              dataKey="meanEnergy"
              stroke="#7de0b2"
              strokeWidth={2}
              dot={false}
              isAnimationActive={false}
            />
          </LineChart>
        </ResponsiveContainer>
      </div>
    </section>
  );
}
