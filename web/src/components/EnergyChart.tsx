import { useMemo } from 'react';
import {
  Area,
  CartesianGrid,
  ComposedChart,
  Legend,
  Line,
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
          <ComposedChart
            data={data}
            margin={{ top: 8, right: 12, left: 8, bottom: 8 }}
          >
            <defs>
              <linearGradient id="localStrokeGradient" x1="0" y1="0" x2="1" y2="0">
                <stop offset="0%" stopColor="#2db9ff" />
                <stop offset="100%" stopColor="#8ed6ff" />
              </linearGradient>
              <linearGradient id="meanStrokeGradient" x1="0" y1="0" x2="1" y2="0">
                <stop offset="0%" stopColor="#59e6b0" />
                <stop offset="100%" stopColor="#b8f7de" />
              </linearGradient>
              <linearGradient id="localAreaGradient" x1="0" y1="0" x2="0" y2="1">
                <stop offset="0%" stopColor="#57b8ff" stopOpacity={0.24} />
                <stop offset="100%" stopColor="#57b8ff" stopOpacity={0.02} />
              </linearGradient>
              <linearGradient id="meanAreaGradient" x1="0" y1="0" x2="0" y2="1">
                <stop offset="0%" stopColor="#7de0b2" stopOpacity={0.18} />
                <stop offset="100%" stopColor="#7de0b2" stopOpacity={0.01} />
              </linearGradient>
            </defs>

            <CartesianGrid strokeDasharray="2 4" stroke="#253749" />
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
              cursor={{ stroke: '#9fb5d8', strokeOpacity: 0.45, strokeWidth: 1 }}
              contentStyle={{
                backgroundColor: '#0f1725',
                border: '1px solid #355070',
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
                stroke="#ffd56a"
                strokeWidth={2}
                strokeDasharray="6 4"
              />
            ) : null}
            {currentPoint ? (
              <>
                <ReferenceDot
                  x={currentPoint.step}
                  y={currentPoint.localEnergy}
                  r={4}
                  fill="#56b5ff"
                  stroke="#d6ebff"
                  strokeWidth={1.5}
                />
                <ReferenceDot
                  x={currentPoint.step}
                  y={currentPoint.meanEnergy}
                  r={4}
                  fill="#7de0b2"
                  stroke="#d8ffec"
                  strokeWidth={1.5}
                />
              </>
            ) : null}
            <Area
              type="monotone"
              dataKey="localEnergy"
              fill="url(#localAreaGradient)"
              stroke="none"
              isAnimationActive={false}
            />
            <Area
              type="monotone"
              dataKey="meanEnergy"
              fill="url(#meanAreaGradient)"
              stroke="none"
              isAnimationActive={false}
            />
            <Line
              type="monotone"
              dataKey="localEnergy"
              stroke="url(#localStrokeGradient)"
              strokeWidth={2}
              dot={false}
              activeDot={{ r: 3, strokeWidth: 0, fill: '#8ed6ff' }}
              isAnimationActive={false}
            />
            <Line
              type="monotone"
              dataKey="meanEnergy"
              stroke="url(#meanStrokeGradient)"
              strokeWidth={2}
              dot={false}
              activeDot={{ r: 3, strokeWidth: 0, fill: '#b8f7de' }}
              isAnimationActive={false}
            />
          </ComposedChart>
        </ResponsiveContainer>
      </div>
    </section>
  );
}
