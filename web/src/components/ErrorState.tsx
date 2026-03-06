interface ErrorStateProps {
  message: string;
  onRetry?: () => void;
}

export function ErrorState({ message, onRetry }: ErrorStateProps) {
  return (
    <div className="error-state">
      <div className="state-card state-card-error">
        <h1 className="state-title">Replay Load Failed</h1>
        <p className="state-text">{message}</p>
        {onRetry ? (
          <button className="action-btn" type="button" onClick={onRetry}>
            Retry
          </button>
        ) : null}
      </div>
    </div>
  );
}
