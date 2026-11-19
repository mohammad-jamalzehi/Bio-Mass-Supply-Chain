import React, { useState } from 'react'
import source from './content/results.json'

const nodes = [
  ['wc', 'Wood chips', 15, 28], ['wp', 'Wood pellets', 15, 103], ['bq', 'Briquettes', 15, 178],
  ['receiving', 'Receiving', 170, 103], ['storage', 'Storage', 330, 178],
  ['existing', 'Existing plant', 485, 28], ['new', 'New plant', 485, 103],
  ['gas', 'Natural gas', 485, 178], ['demand', 'Campus heat', 650, 103],
] as const
const edges = [
  { from: 'wc', to: 'receiving', path: 'M141 52 H151 Q156 52 156 60 V119 Q156 127 170 127' },
  { from: 'wp', to: 'receiving', path: 'M141 127 H170' },
  { from: 'bq', to: 'receiving', path: 'M141 202 H151 Q156 202 156 194 V135 Q156 127 170 127' },
  { from: 'receiving', to: 'storage', path: 'M296 135 H310 Q315 135 315 145 V192 Q315 202 330 202' },
  { from: 'receiving', to: 'existing', path: 'M296 117 H326 Q341 117 341 102 V66 Q341 52 358 52 H485' },
  { from: 'receiving', to: 'new', path: 'M296 127 H485' },
  { from: 'storage', to: 'existing', path: 'M377 178 V72 Q377 63 389 63 H485' },
  { from: 'storage', to: 'new', path: 'M455 194 H460 Q466 194 466 160 V146 Q466 138 485 138' },
  { from: 'existing', to: 'demand', path: 'M611 52 H625 Q632 52 632 66 V115 Q632 127 650 127', heat: true },
  { from: 'new', to: 'demand', path: 'M611 127 H650', heat: true },
  { from: 'gas', to: 'demand', path: 'M611 202 H625 Q632 202 632 190 V140 Q632 127 650 127', heat: true },
]
const caseTitles: Record<string, string> = { base: 'Base storage', 'tight-silo': 'Tighter pellet silo', 'no-storage': 'No day-one storage' }
const fixtures = source.fixtures.filter((item) => item.phase === 3)
const value = (n: number) => n.toLocaleString('en-US', { maximumFractionDigits: 3 })

export function Network({ preview = false }: { preview?: boolean }) {
  const [index, setIndex] = useState(0)
  const [dayIndex, setDay] = useState(0)
  const [selected, setSelected] = useState('storage')
  const fixture = fixtures[index]
  const day = fixture.days[dayIndex]
  if (!('burns' in day)) throw new Error('Expected a multi-plant fixture')
  const inventory = day.inventory_by_fuel
  const burns = day.burns
  const heat = day.biomass_gj_by_plant_and_fuel
  const fluids = (kind: 'fresh_m3' | 'stored_m3', plant?: string) => Object.entries(burns ?? {}).filter(([p]) => !plant || p === plant)
    .flatMap(([, fuels]) => Object.values(fuels).map((flow) => flow[kind])).reduce((sum, n) => sum + n, 0)
  const delivered = Object.values(inventory ?? {}).reduce((sum, item) => sum + item.delivered, 0)
  const closing = Object.values(inventory ?? {}).reduce((sum, item) => sum + item.closing, 0)
  const useful = (plant: string) => Object.values(heat?.[plant as keyof typeof heat] ?? {}).reduce((sum, n) => sum + n, 0)
  const gas = day.gas_gj ?? 0
  const demand = useful('existing') + useful('new') + gas
  const storageIn = Math.max(0, delivered - fluids('fresh_m3'))
  function volume(from: string, to: string) {
    if (from === 'gas') return gas
    if (to === 'demand') return useful(from)
    if (to === 'receiving') return inventory?.[from as keyof typeof inventory]?.delivered ?? 0
    if (to === 'storage') return storageIn
    return fluids(from === 'storage' ? 'stored_m3' : 'fresh_m3', to)
  }
  function detail() {
    const stock = inventory?.[selected as keyof typeof inventory]
    if (stock) return `Delivered ${value(stock.delivered)} m³ · closing stock ${value(stock.closing)} m³`
    if (selected === 'storage') return `Fresh receipts kept ${value(storageIn)} m³ · withdrawn ${value(fluids('stored_m3'))} m³ · closing stock ${value(closing)} m³`
    if (selected === 'receiving') return `${value(delivered)} m³ delivered · ${value(fluids('fresh_m3'))} m³ sent directly to plants`
    if (selected === 'gas') return `${value(gas)} GJ useful gas heat`
    if (selected === 'demand') return `${value(demand)} GJ met by biomass and gas`
    return `${value(fluids('fresh_m3', selected))} m³ fresh · ${value(fluids('stored_m3', selected))} m³ stored · ${value(useful(selected))} GJ useful output`
  }
  const marker = preview ? 'hero-arrow' : 'network-arrow'
  return <div className={preview ? 'preview-network' : 'network-card'}>
    {!preview && <div className="graph-controls"><span className="pill">Synthetic solver data</span>
      <label>Scenario <select value={index} onChange={(event) => { setIndex(Number(event.target.value)); setDay(0) }}>
        {fixtures.map((item, i) => <option key={item.scenario} value={i}>{caseTitles[item.scenario]}</option>)}
      </select></label>
      <label>Day <select value={dayIndex} onChange={(event) => setDay(Number(event.target.value))}>
        {fixture.days.map((_, i) => <option key={i} value={i}>{i + 1}</option>)}
      </select></label></div>}
    <div className="graph-scroll" tabIndex={preview ? undefined : 0} role="region" aria-label={preview ? 'System topology' : 'Scrollable optimization network'}>
      <svg className="diagram" viewBox="0 0 790 242" role={preview ? 'img' : undefined} aria-label="Three feedstocks enter receiving. Fresh and stored biomass supply both plants. Biomass heat and gas meet demand.">
        <defs><marker id={marker} viewBox="0 0 10 10" refX="8" refY="5" markerWidth="6" markerHeight="6" orient="auto"><path d="M0 0 L10 5 L0 10 Z" fill="context-stroke" /></marker></defs>
        {edges.map((edge) => {
          const actual = volume(edge.from, edge.to)
          const active = !preview && actual > 1e-8
          return <path key={edge.from + edge.to} d={edge.path} markerEnd={`url(#${marker})`}
            className={`flow ${edge.heat ? 'heat' : ''} ${active ? 'has-flow' : ''}`}><title>{preview ? (edge.heat ? 'Useful heat route' : 'Material route') : `${value(actual)} ${edge.heat ? 'GJ useful' : 'm³'}`}</title></path>
        })}
        {nodes.map(([id, title, x, y]) => <g key={id} className={`graph-node ${!preview && selected === id ? 'node-selected' : ''}`}
          role={preview ? undefined : 'button'} tabIndex={preview ? undefined : 0} aria-label={preview ? undefined : title}
          aria-pressed={preview ? undefined : selected === id} onClick={preview ? undefined : () => setSelected(id)}
          onKeyDown={preview ? undefined : (event) => { if (event.key === 'Enter' || event.key === ' ') { event.preventDefault(); setSelected(id) } }}>
          <rect x={x} y={y} width="126" height="48" rx="8" /><text x={x + 63} y={y + 29} textAnchor="middle">{title}</text></g>)}
      </svg>
    </div>
    {!preview && <><div className="graph-legend"><span>Biomass volume (m³)</span><span>Useful heat (GJ)</span><span>Bold paths carry positive flow</span></div>
      <div className="graph-detail" aria-live="polite"><strong>{nodes.find(([id]) => id === selected)?.[1]}</strong><span>{detail()}</span></div>
      <p className="network-note">Select a node or inspect an edge for a quantity. Material flows and energy flows use independent units. Paths indicate routing, with line weight distinguishing positive flow from zero.</p></>}
  </div>
}
