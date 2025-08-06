# PharmaPy CI/CD Configuration

This repository uses a consolidated GitHub Actions workflow to ensure efficient testing and reduce CI resource consumption.

## Workflow Structure

### `consolidated-ci.yml`
A single, comprehensive workflow that replaces multiple separate workflows, reducing the total number of CI jobs from 99+ to approximately 12-15.

#### Jobs:

1. **Core Tests** (6 jobs)
   - Ubuntu + Windows × Python 3.9, 3.11
   - macOS × Python 3.11
   - Ubuntu × Python 3.12
   - Runs unit and integration tests
   - Tests installation methods

2. **Assimulo Integration** (1 job)
   - Ubuntu 20.04 × Python 3.9
   - Tests scientific computing integration
   - Specialized for numerical solver compatibility

3. **Installation Scripts** (2 jobs, conditional)
   - Ubuntu + Windows
   - Only runs on main branches
   - Tests installation scripts

4. **Code Quality & Documentation** (1 job, conditional)
   - Only runs on main/master branches
   - Linting, coverage, documentation build

5. **Build & Package** (1 job, conditional)
   - Only runs on main/master branches or tags
   - Creates distributable packages

## Triggers

- **Push**: main, master, develop, test-devel branches
- **Pull Request**: main, master, develop branches  
- **Schedule**: Weekly on Sundays at 3 AM UTC
- **Manual**: Via workflow_dispatch

## Resource Optimization

- **Reduced OS Matrix**: Focus on Ubuntu/Windows, selective macOS
- **Strategic Python Versions**: LTS (3.9) + Current (3.11) + Latest (3.12)
- **Conditional Jobs**: Quality checks and builds only on main branches
- **Smart Matrix Includes**: Additional combinations only where needed

This configuration provides comprehensive testing coverage while using ~85% fewer CI resources than the previous setup.
