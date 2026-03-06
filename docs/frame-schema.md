# Message type examples

## Init

```json
{
  "type": "init",
  "numParticles": 54,
  "boxLength": 10.0,
  "warmupSteps": 5000,
  "measureSteps": 50000,
  "stepSize": 0.35,
  "blockSize": 100
}
```

## Frame

```json
{
  "type": "frame",
  "step": 120,
  "accepted": 67,
  "proposed": 120,
  "acceptanceRate": 0.5583,
  "localEnergy": -1.234567,
  "meanEnergy": -1.198234,
  "standardError": 0.0123,
  "positions": [1.2, 3.4, 5.6, 2.1, 8.0, 0.2]
}
```

## Type

```json
{
  "type": "done",
  "finalAcceptanceRate": 0.5086,
  "finalMeanEnergy": -1.20311,
  "finalStandardError": 0.00482
}
```
