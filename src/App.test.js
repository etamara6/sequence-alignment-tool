import { render, screen } from '@testing-library/react';
import App from './App';

test('renders the main app title', () => {
  render(<App />);
  const titleElement = screen.getByText(/Sequence Alignment Studio/i);
  expect(titleElement).toBeInTheDocument();
});
