# Design System
## High-Throughput Molecular Docking Pipeline — Web Dashboard

**Version:** 1.0  
**Last Updated:** 2026-11-07  
**Scope:** Optional web dashboard (Phase 6+). The CLI pipeline has no visual design requirements.

> This document defines what the design IS, not just what it avoids. Refer to `rules.md` Part D for the full banned-patterns list.

---

## 1. Design Philosophy

**Scientific clarity over marketing aesthetics.**  
This is a research tool. Every visual decision should make data easier to interpret, not make the page look impressive. If a design element cannot be justified by usability or data clarity, remove it.

**Hierarchy through typography and whitespace, not colour.**  
Use font weight, size, and spacing to establish hierarchy. Colour is used sparingly as a semantic signal (e.g., affinity score heat, error states), never for decoration.

**Dense but not cramped.**  
Scientific users expect data-dense interfaces. Tables and charts should show as much meaningful data as possible without sacrificing readability. Use fixed-width columns for numeric data.

---

## 2. Colour Palette

### Light Mode (default)

| Token | Hex | Use |
|---|---|---|
| `--bg-base` | `#FAFAF9` | Page background (warm off-white, not pure white) |
| `--bg-surface` | `#FFFFFF` | Card and panel backgrounds |
| `--bg-subtle` | `#F4F4F3` | Table row alternating, input backgrounds |
| `--border` | `#E4E4E0` | All borders |
| `--text-primary` | `#1A1A18` | Body text, headings |
| `--text-secondary` | `#6B6B63` | Labels, metadata, secondary text |
| `--text-muted` | `#9C9C94` | Placeholder text, disabled states |
| `--accent` | `#2563EB` | Primary action colour (single blue, no gradient) |
| `--accent-hover` | `#1D4ED8` | Accent hover state (darker, not transparent) |
| `--success` | `#16A34A` | Success states, positive affinity |
| `--warning` | `#D97706` | Warnings, partial failures |
| `--error` | `#DC2626` | Errors, failed docking |

### Dark Mode

| Token | Hex | Use |
|---|---|---|
| `--bg-base` | `#111110` | Page background |
| `--bg-surface` | `#1C1C1A` | Cards, panels |
| `--bg-subtle` | `#232320` | Table rows, input backgrounds |
| `--border` | `#2E2E2B` | All borders |
| `--text-primary` | `#EDEDE8` | Body text |
| `--text-secondary` | `#9C9C94` | Secondary text |
| `--text-muted` | `#6B6B63` | Muted text |
| `--accent` | `#3B82F6` | Primary action |
| `--accent-hover` | `#60A5FA` | Accent hover |

**Dark mode contrast requirement:** Text on background must meet WCAG AA (minimum 4.5:1 for body, 3:1 for large text). No low-contrast combinations.

---

## 3. Typography

### Font Stack

```css
--font-sans: 'DM Sans', -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif;
--font-mono: 'JetBrains Mono', 'Fira Code', 'Cascadia Code', monospace;
```

**DM Sans** for all UI text. Geometric but warm. Legible at small sizes. Scientific without being cold.  
**JetBrains Mono** for all code, SMILES strings, PDB IDs, residue names, file paths.

No serif fonts. No display fonts. No mixing of two sans-serif families.

### Type Scale

| Token | Size | Weight | Use |
|---|---|---|---|
| `--text-xs` | 11px | 400 | Table footnotes, metadata |
| `--text-sm` | 13px | 400 | Table cells, secondary labels |
| `--text-base` | 15px | 400 | Body text, default |
| `--text-lg` | 17px | 500 | Section labels, card titles |
| `--text-xl` | 20px | 600 | Panel headings |
| `--text-2xl` | 24px | 700 | Page titles |

Numeric data in tables: use `font-variant-numeric: tabular-nums` so digits align in columns.

---

## 4. Spacing System

Base unit: **4px**. All spacing is a multiple of 4.

| Token | Value | Use |
|---|---|---|
| `--space-1` | 4px | Tight intra-element gaps |
| `--space-2` | 8px | Icon-to-label, tight padding |
| `--space-3` | 12px | Default component padding |
| `--space-4` | 16px | Card padding, section gaps |
| `--space-6` | 24px | Between cards |
| `--space-8` | 32px | Section headers |
| `--space-12` | 48px | Between major page sections |
| `--space-16` | 64px | Page-level vertical rhythm |

No `em` or `rem` mixing. Use `px` for components, `rem` for type only.

---

## 5. Component Patterns

### Cards

```
Background: var(--bg-surface)
Border: 1px solid var(--border)
Border-radius: 8px
Padding: var(--space-4)
Shadow: 0 1px 3px rgba(0,0,0,0.06)
```

No coloured borders. No glassmorphism (backdrop-filter: blur). The card surface IS its own layer — don't fake depth with blur.

### Tables (primary UI for docking results)

- Fixed-width columns for numeric data (affinity, RMSD)
- Monospace font for SMILES, residue names, PDB IDs
- Row hover: `--bg-subtle` background change — no colour change
- Sortable columns: show sort direction with a simple caret, no animation
- Paginate at 50 rows; show total count

### Buttons

| Type | Background | Text | Hover |
|---|---|---|---|
| Primary | `--accent` | white | `--accent-hover` (solid) |
| Secondary | `--bg-subtle` | `--text-primary` | `--border` background |
| Destructive | `#FEF2F2` | `--error` | `#FEE2E2` background |

Hover state changes background to a different solid colour. No opacity fades.

### Status Badges

| Status | Background | Text |
|---|---|---|
| `success` | `#F0FDF4` | `#166534` |
| `prep_failed` | `#FFFBEB` | `#92400E` |
| `docking_failed` | `#FEF2F2` | `#991B1B` |

Badges use a 4px border-radius and 6px/12px padding. No rounded-full "pill" badges for status.

### Charts

**Affinity distribution:** Histogram using Recharts `BarChart`. Single colour (`--accent`). No gradients on bars.  
**ROC curve:** Line chart. Single line, `--accent` colour. Reference diagonal in `--border` colour. No fill under curve.  
**Top-5 hits table:** Plain `<table>` styled with the table spec above. Not a card grid.

---

## 6. Layout

### Dashboard Layout

```
┌─────────────────────────────────────┐
│ Header: logo + run status           │
├──────────────┬──────────────────────┤
│ Sidebar      │ Main content area    │
│ (nav only)   │                      │
│              │                      │
└──────────────┴──────────────────────┘
```

Sidebar: 240px fixed width. Navigation only. No collapsible panels in the sidebar.  
Main area: flexible. Max-width 1200px. Centered on wide screens.

### Page Sections

Sections are separated by whitespace (`--space-12`), not dividers or horizontal rules. Use heading hierarchy to signal section boundaries.

---

## 7. Iconography

Use **Phosphor Icons** (not Lucide). Import only the specific icons you need. Icons are used exclusively for:

- Navigation (sidebar items)
- Action buttons (copy, download, sort)
- Status indicators (success/warning/error)

Icons must never replace text labels. Always pair icon + label for primary navigation. Standalone icons only in tight spaces with a tooltip.

Size: 16px for inline icons, 20px for standalone action icons.

---

## 8. Accessibility Requirements

- All text meets WCAG AA contrast (4.5:1 body, 3:1 large)
- All interactive elements reachable via keyboard
- Focus ring: `2px solid --accent` with `2px offset` — never `outline: none` without a custom replacement
- Tables use `<thead>`, `<th scope="col">`, `<caption>` correctly
- Charts have aria-label descriptions
- No information conveyed by colour alone — always pair with text or shape
