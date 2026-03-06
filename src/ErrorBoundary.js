import { Component } from "react";

/**
 * Catches unexpected runtime errors inside the alignment UI
 * and shows a friendly fallback instead of a blank/broken page.
 */
class ErrorBoundary extends Component {
  constructor(props) {
    super(props);
    this.state = { hasError: false, message: "" };
  }

  static getDerivedStateFromError(error) {
    return { hasError: true, message: error.message };
  }

  componentDidCatch(error, info) {
    // In production you'd send this to a logging service (e.g. Sentry)
    console.error("Alignment error:", error, info);
  }

  render() {
    if (this.state.hasError) {
      return (
        <div style={{
          padding: "2rem",
          textAlign: "center",
          color: "#880e4f",
          fontFamily: "monospace"
        }}>
          <h2>⚠️ Something went wrong</h2>
          <p>{this.state.message}</p>
          <button onClick={() => this.setState({ hasError: false, message: "" })}>
            Try again
          </button>
        </div>
      );
    }
    return this.props.children;
  }
}

export default ErrorBoundary;