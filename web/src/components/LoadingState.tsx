export function LoadingState() {
  return (
    <div className="loading-state">
      <div className="state-card">
        <h1 className="state-title">Loading Replay</h1>
        <p className="state-text">Parsing /sample-run.jsonl and preparing frames...</p>
      </div>
    </div>
  );
}
