.PHONY: help init install-uv install dev conda-setup test test-fast test-unit test-int lint lint-fix format format-check type-check check ci-check clean generate-testdata lock sync

help:  ## Show this help message
	@echo "Usage: make [target]"
	@echo ""
	@echo "Targets:"
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | \
		awk 'BEGIN {FS = ":.*?## "}; {printf "  \033[36m%-20s\033[0m %s\n", $$1, $$2}'

# ==================== INSTALLATION ====================

init: install-uv dev  ## Initialize complete development environment

install-uv:  ## Install uv package manager (if not present)
	@command -v uv >/dev/null 2>&1 || { \
		echo "Installing uv..."; \
		curl -LsSf https://astral.sh/uv/install.sh | sh; \
	}
	@echo "uv installed"

install:  ## Install package in production mode
	uv pip install .

dev:  ## Install package with development dependencies
	uv pip install -e ".[dev]"
	@echo "Development environment ready"

conda-setup:  ## Create conda environment for bioinformatics tools
	@command -v mamba >/dev/null 2>&1 || command -v conda >/dev/null 2>&1 || { \
		echo "ERROR: conda/mamba not found."; \
		exit 1; \
	}
	@if command -v mamba >/dev/null 2>&1; then \
		mamba env create -f conda/environment.yml --force; \
	else \
		conda env create -f conda/environment.yml --force; \
	fi
	@echo "Conda environment created: open-pacmuci-tools"

# ==================== TESTING ====================

test:  ## Run all tests with coverage
	uv run pytest

test-fast:  ## Run unit tests only, no coverage (fast)
	uv run pytest tests/unit/ --no-cov -x

test-unit:  ## Run unit tests with coverage
	uv run pytest tests/unit/

test-int:  ## Run integration tests only
	uv run pytest tests/integration/ -m "integration or e2e"

# ==================== CODE QUALITY ====================

lint:  ## Run ruff linter
	uv run ruff check src/open_pacmuci/ tests/

lint-fix:  ## Run ruff linter and auto-fix issues
	uv run ruff check --fix src/open_pacmuci/ tests/

format:  ## Format code with ruff
	uv run ruff format src/open_pacmuci/ tests/

format-check:  ## Check if code is formatted (no changes)
	uv run ruff format --check src/open_pacmuci/ tests/

type-check:  ## Run mypy type checker
	uv run mypy src/open_pacmuci/

check: lint format-check type-check test  ## Run all quality checks

ci-check:  ## Run EXACT same checks as GitHub Actions CI
	@echo "Running CI checks locally..."
	@echo ""
	@echo "=== Code Quality Checks ==="
	@echo "1. Ruff linter..."
	uv run ruff check src/open_pacmuci/ tests/
	@echo "2. Ruff formatter check..."
	uv run ruff format --check src/open_pacmuci/ tests/
	@echo "3. Mypy type checker..."
	uv run mypy src/open_pacmuci/
	@echo ""
	@echo "=== Test Suite ==="
	@echo "4. Running pytest with coverage (>=80%)..."
	uv run pytest tests/unit/ --cov=open_pacmuci --cov-report=term-missing --cov-fail-under=80
	@echo ""
	@echo "All CI checks passed!"

# ==================== TEST DATA ====================

generate-testdata:  ## Generate test data via MucOneUp (requires conda env + MucOneUp)
	uv run python scripts/generate_testdata.py

# ==================== CLEANUP ====================

clean:  ## Remove build artifacts and caches
	rm -rf build/ dist/ *.egg-info
	find . -type d -name __pycache__ -exec rm -rf {} + 2>/dev/null || true
	find . -type f -name "*.pyc" -delete
	find . -type d -name ".pytest_cache" -exec rm -rf {} + 2>/dev/null || true
	find . -type d -name ".mypy_cache" -exec rm -rf {} + 2>/dev/null || true
	find . -type d -name ".ruff_cache" -exec rm -rf {} + 2>/dev/null || true
	find . -type d -name "htmlcov" -exec rm -rf {} + 2>/dev/null || true
	rm -f .coverage coverage.xml
	@echo "Cleaned build artifacts"

# ==================== UTILITIES ====================

lock:  ## Update uv.lock file
	uv lock

sync:  ## Sync environment with uv.lock
	uv sync

.DEFAULT_GOAL := help
