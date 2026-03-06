import { useMemo } from 'react';
import {
  CartesianGrid,
  Legend,
  Line,
  LineChart,
  ReferenceDot,
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

  const currentPoint = frames[currentFrameIndex];
  const currentStep = currentPoint?.step;

  return (
    <section className="panel-card chart-panel">
      <h2 className="panel-title">Energy History</h2>
      <div className="chart-shell">
        <ResponsiveContainer width="100%" height="100%">
          <LineChart data={data} margin={{ top: 8, right: 12, left: 8, bottom: 8 }}>
            <CartesianGrid strokeDasharray="3 3" stroke="#1b222c" />
            <XAxis
              dataKey="step"
              tick={{ fill: '#858c99', fontSize: 11 }}
              axisLine={{ stroke: '#252d39' }}
              tickLine={{ stroke: '#252d39' }}
            />
            <YAxis
              tick={{ fill: '#858c99', fontSize: 11 }}
              axisLine={{ stroke: '#252d39' }}
              tickLine={{ stroke: '#252d39' }}
              domain={['auto', 'auto']}
            />
            <Tooltip
              cursor={{ stroke: '#5d6674', strokeOpacity: 0.6, strokeWidth: 1 }}
              contentStyle={{
                backgroundColor: '#0b0d11',
                border: '1px solid #252d39',
                borderRadius: 6,
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
              wrapperStyle={{ color: '#bfc6d1', fontSize: 11 }}
            />
            {Number.isFinite(currentStep) ? (
              <ReferenceLine
                x={currentStep}
                stroke="#b8a368"
                strokeWidth={1.5}
                strokeDasharray="4 4"
              />
            ) : null}
            {currentPoint ? (
              <>
                <ReferenceDot
                  x={currentPoint.step}
                  y={currentPoint.localEnergy}
                  r={3}
                  fill="#57b8ff"
                  stroke="none"
                />
                <ReferenceDot
                  x={currentPoint.step}
                  y={currentPoint.meanEnergy}
                  r={3}
                  fill="#74d1a7"
                  stroke="none"
                />
              </>
            ) : null}
            <Line
              type="monotone"
              dataKey="localEnergy"
              stroke="#57b8ff"
              strokeWidth={1.9}
              dot={false}
              isAnimationActive={false}
            />
            <Line
              type="monotone"
              dataKey="meanEnergy"
              stroke="#74d1a7"
              strokeWidth={1.9}
              dot={false}
              isAnimationActive={false}
            />
          </LineChart>
        </ResponsiveContainer>
      </div>
    </section>
  );
}
