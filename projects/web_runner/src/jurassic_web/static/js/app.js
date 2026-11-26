/* JURASSIC Web Runner · no front-end framework required, no accounts or tracking. */
(() => {
  'use strict';
  const $ = (id) => document.getElementById(id);
  const $$ = (sel) => Array.from(document.querySelectorAll(sel));
  const csrf = document.querySelector('meta[name="csrf-token"]').content;
  const state = { catalog: {}, dataset: 'airs', geometry: 'nadir', gases: ['CO2'], channels: [], runs: [],
                  active: null, view: 'plots', metric: 'radiance', section: 'simulation', poller: null, output: 'bbt',
                  plotMode: 'image', plotRequest: 0 };
  const escapeHtml = (x) => String(x).replace(/[&<>"']/g, (v) => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[v]));
  const text = (id, value) => { if ($(id)) $(id).textContent = String(value); };
  let toastTimer;
  const toast = (message, bad=false) => { const node=$('toast'); node.textContent=message; node.classList.toggle('error',bad); node.classList.add('visible');
    clearTimeout(toastTimer); toastTimer=setTimeout(()=>node.classList.remove('visible'),5500); };
  async function api(url, method='GET', data) {
    const opt = {method,headers:{'X-CSRF-Token':csrf}};
    if (data instanceof FormData) opt.body=data;
    else if (data !== undefined) { opt.body=JSON.stringify(data); opt.headers['Content-Type']='application/json'; }
    let response;
    try { response=await fetch(url,opt); } catch(e) {throw new Error('Cannot connect to the web runner server.');}
    const content = await response.text();
    let result;
    try {result=JSON.parse(content);} catch(e) { throw new Error('The server returned an unexpected response.'); }
    if (!response.ok) throw new Error(result.error || `Server error ${response.status}`);
    return result;
  }
  function theme(name) {
    document.documentElement.dataset.theme=name;
    localStorage.setItem('jurassic-theme',name);
    $('light-btn').classList.toggle('selected',name==='light');
    $('dark-btn').classList.toggle('selected',name==='dark');
    if (state.section === 'simulation') renderMetric();
    if (state.section === 'compare' && $('compare-plot').dataset.ready) doCompare().catch(()=>{});
  }
  function section(name) {
    if (!$(name+'-section')) return;
    state.section=name;
    $$('.section').forEach(el => {const show=el.id===name+'-section';el.hidden=!show;el.classList.toggle('active-section',show);});
    $$('.nav-item').forEach(el=>el.classList.toggle('active',el.dataset.section===name && (name!=='datasets'||el.textContent.trim()==='Datasets')));
    if(name==='history') renderHistory();
    if(name==='datasets') renderCatalog();
    if(name==='compare') renderComparePicker();
    if(name==='settings') refreshSystem();
    if(name==='configurations') previewConfig().then(()=>text('config-studio-text',$('config-text').textContent)).catch(()=>{});
    if(name==='simulation') setTimeout(()=>renderMetric(),30);
  }
  function flow(name) {
    if(name==='results') {section('simulation'); setView('plots'); $('results-title').scrollIntoView({behavior:'smooth',block:'nearest'});}
    else if(name==='run') {section('simulation'); $('run-button').focus();}
    else {section('simulation'); $('setup-title').scrollIntoView({behavior:'smooth',block:'nearest'});}
    $$('.flow-step').forEach(el=>el.classList.toggle('active',el.dataset.flow===name));
  }
  function setView(name) {
    state.view=name;
    $$('#content-tabs button').forEach(el=>el.classList.toggle('selected',el.dataset.view===name));
    ['plots','table','atmosphere','config','logs'].forEach(id=>$(id+'-view').hidden=id!==name);
    if(name==='plots') renderMetric();
    if(name==='table') renderTable();
    if(name==='atmosphere') renderAtmosphere();
    if(name==='logs') renderLogs();
    if(name==='config') renderConfig();
  }
  function currentDataset() {return state.catalog[state.dataset] || null;}
  function allChannels(ds) {return ds?Object.keys(ds.channels).sort((a,b)=>Number(a)-Number(b)):[];}
  function allGases(ds) {return ds?[...new Set(allChannels(ds).flatMap(c=>Object.keys(ds.channels[c].gases)))].sort():[];}
  function defaultOutput(geo) {return geo==='limb'?'radiance':'bbt';}
  function setGeometry(geo) {
    state.geometry=geo;
    state.dataset = ({nadir:'airs',limb:'boxcar',zenith:'boxcar'})[geo] in state.catalog
                  ? ({nadir:'airs',limb:'boxcar',zenith:'boxcar'})[geo]
                  : Object.keys(state.catalog)[0] || '';
    state.output=defaultOutput(geo);
    $$('input[name="output"]').forEach(e=>e.checked=e.value===state.output);
    resetSelections(); paintSetup();
  }
  function resetSelections() {
    const ds=currentDataset();
    state.channels=allChannels(ds);
    const gasList=allGases(ds);
    state.gases=state.dataset==='airs'?gasList.filter(g=>g==='CO2'):gasList.filter(g=>['CO2','H2O','O3'].includes(g));
    if(!state.gases.length && gasList.length) state.gases=[gasList[0]];
  }
  function paintSetup() {
    $$('#geometry-options button').forEach(el=>el.classList.toggle('selected',el.dataset.geometry===state.geometry));
    const select=$('dataset-select');
    select.replaceChildren();
    Object.entries(state.catalog).forEach(([key,ds])=>{
      const opt=new Option(`${ds.label} (${allChannels(ds).join(' / ')} cm⁻¹)`,key);
      select.add(opt);
    });
    select.value=state.dataset;
    const ds=currentDataset();
    if(!ds){ $('gases-list').textContent='No datasets found';$('channels-list').textContent='No channels'; return; }
    $('channels-list').innerHTML=allChannels(ds).map(c=>`<label><input type="checkbox" data-channel="${escapeHtml(c)}" ${state.channels.includes(c)?'checked':''}><span>${escapeHtml(c)} cm⁻¹</span></label>`).join('');
    $('gases-list').innerHTML=allGases(ds).map(g=>{
      const available=state.channels.some(ch=>ds.channels[ch].gases[g]?.data);
      return `<label class="${available?'':'disabled'}"><input type="checkbox" data-gas="${escapeHtml(g)}" ${state.gases.includes(g)?'checked':''} ${available?'':'disabled'}><span>${escapeHtml(g)}${available?'':' · no data'}</span></label>`;
    }).join('');
    $$('#channels-list input').forEach(el=>el.addEventListener('change',()=>{
      state.channels=$$('#channels-list input:checked').map(x=>x.dataset.channel);paintSetup();
    }));
    $$('#gases-list input').forEach(el=>el.addEventListener('change',()=>{
      state.gases=$$('#gases-list input:checked').map(x=>x.dataset.gas);refreshCompatibility();
    }));
    refreshCompatibility();
  }
  function refreshCompatibility(){
    const ds=currentDataset(); if(!ds)return;
    const warnings=[];
    if(!state.gases.length||!state.channels.length)warnings.push('Select at least one gas and one channel.');
    for(const ch of state.channels) {
      const item=ds.channels[ch];if(!item)continue;
      for(const gas of state.gases) if(!item.gases[gas]?.data) warnings.push(`${gas} at ${ch} cm⁻¹ has no numerical lookup-table data.`);
    }
    if(state.gases.some(g=>!['CO2','H2O','O3','F11','CCl4'].includes(g)))warnings.push('Custom gas selected: upload a matching atmosphere in Advanced Mode.');
    const note=$('compatibility-note');
    note.classList.toggle('warning',warnings.length>0);
    note.textContent=warnings.length?warnings.slice(0,3).join('  ·  '):'Dataset choices have matching numeric lookup-table files. Scientific suitability is not guaranteed.';
  }
  function parameters() {
    return {geometry:state.geometry,dataset:state.dataset,gases:state.gases,channels:state.channels,
      output:document.querySelector('input[name="output"]:checked').value,
      rayds:parseFloat($('rayds').value),raydz:parseFloat($('raydz').value),
      refrac:$('refrac').checked,kernel:$('kernel').checked};
  }
  async function previewConfig() {
    const data=await api('/api/config/preview','POST',parameters());
    $('config-text').textContent=data.text;
    $('config-studio-text').textContent=data.text;
    return data;
  }
  function saveFile(content, name, type='text/plain') {
    const url=URL.createObjectURL(new Blob([content],{type}));
    const a=document.createElement('a');a.href=url;a.download=name;document.body.appendChild(a);a.click();a.remove();
    setTimeout(()=>URL.revokeObjectURL(url),1000);
  }
  async function runSimulation() {
    const runButton=$('run-button');
    runButton.disabled=true;runButton.querySelector('span').textContent='Preparing simulation…';
    try {
      const form=new FormData();form.append('payload',JSON.stringify(parameters()));
      const file=$('atmosphere-file').files[0];if(file)form.append('atmosphere',file);
      const item=await api('/api/runs','POST',form);
      state.active=item.id;
      text('run-status-pill','Queued');$('run-status-pill').className='status-pill running';
      text('run-duration','Running…');
      toast('Simulation queued. Scientific commands will run in an isolated folder.');
      flow('results');
      await refreshHistory();
      clearInterval(state.poller);
      state.poller=setInterval(()=>pollRun(item.id),1700);
      await pollRun(item.id);
    }catch(e){toast(e.message,true);}finally{runButton.disabled=false;runButton.querySelector('span').textContent='Run Simulation';}
  }
  async function pollRun(id){
    try{
      const r=await api('/api/runs/'+id);
      if(state.active!==id){clearInterval(state.poller);return;}
      updateRunSummary(r);
      if(['completed','failed','interrupted'].includes(r.status)){
        clearInterval(state.poller);state.poller=null;
        await refreshHistory();
        if(r.status==='completed'){toast('Simulation completed. Results are ready.');setView('plots');}
        else toast(r.error||'Simulation stopped. See diagnostics.',true);
      }
    }catch(e){clearInterval(state.poller);toast(e.message,true);}
  }
  function updateRunSummary(r) {
    const status=$('run-status-pill');status.textContent=r.status==='running'?r.stage:r.status;
    status.className='status-pill '+(r.status==='completed'?'success':['queued','running'].includes(r.status)?'running':r.status==='failed'?'error':'');
    text('run-duration',r.duration?`${r.duration.toFixed(1)} s`:'In progress');
    text('run-date',new Date(r.created_at).toLocaleString());
    const p=r.parameters;
    $('run-summary').innerHTML=[['Geometry',p.geometry],['Dataset',p.dataset],['Gases',p.gases.join(', ')],['Channels',p.channels.join(', ')+' cm⁻¹']]
      .map(([a,b])=>`<div><dt>${escapeHtml(a)}</dt><dd>${escapeHtml(b)}</dd></div>`).join('');
    const ready=r.status==='completed';
    $('downloads').innerHTML = (ready?[
      ['rad.tab','rad.tab'],['atm.tab','atm.tab'],['obs.tab','obs.tab'],['config.ctl','config.ctl'],
      ['Plots (ZIP)','plots.zip'],['All Results (ZIP)','all_results.zip'],...(p.kernel?[['kernel.tab','kernel.tab']]:[])
    ].map(([caption,name])=>`<a class="btn small ${name==='all_results.zip'?'primary-mini':''}" href="/api/runs/${escapeHtml(r.id)}/download/${encodeURIComponent(name)}">${escapeHtml(caption)}</a>`):['rad.tab','atm.tab','obs.tab','config.ctl','Plots (ZIP)','All Results (ZIP)'].map(x=>`<button class="btn small" disabled>${x}</button>`)).join('');
    text('chart-note',ready?`Output from run ${r.id} · ${r.warnings} warning(s) · JURASSIC exit codes recorded in commands.json`:
      r.error||`${r.stage} · ${r.status}`);
  }
  async function openRun(id) {
    state.active=id;clearInterval(state.poller);
    const r=await api('/api/runs/'+id);updateRunSummary(r);flow('results');
    if(r.status==='completed')setView('plots');
    else if(['running','queued'].includes(r.status)) state.poller=setInterval(()=>pollRun(id),1700);
  }
  async function refreshHistory() {state.runs=await api('/api/runs');if(state.section==='history')renderHistory();if(state.section==='compare')renderComparePicker();}
  function runHtml(r) {
    const info=r.parameters;
    return `<div class="history-row"><span class="status-indicator ${escapeHtml(r.status)}"></span><div><strong>${escapeHtml(info.geometry.toUpperCase())} · ${escapeHtml(info.gases.join(' + '))}</strong>
      <small>${escapeHtml(new Date(r.created_at).toLocaleString())} · ${escapeHtml(r.status)} · ${escapeHtml(r.id)}</small></div>
      <button data-open="${escapeHtml(r.id)}" type="button">Open →</button></div>`;
  }
  function attachOpen(){ $$('[data-open]').forEach(el=>el.onclick=()=>openRun(el.dataset.open).catch(e=>toast(e.message,true))); }
  function renderHistory(){ $('run-history').innerHTML=state.runs.length?state.runs.map(runHtml).join(''):'<div class="empty-state">No runs yet. Start your first simulation from Setup.</div>';attachOpen(); }
  function renderComparePicker() {
    const finished=state.runs.filter(r=>r.status==='completed');
    $('compare-picker').innerHTML=finished.length?finished.map(r=>
      `<label class="history-row"><input type="checkbox" class="compare-choice" value="${r.id}"><span><strong>${escapeHtml(r.geometry)} · ${escapeHtml(r.parameters.gases.join('+'))}</strong><small>${escapeHtml(r.id)} · ${escapeHtml(r.created_at)}</small></span><span>✓</span></label>`).join(''):
      '<div class="empty-state">Complete at least two runs to compare results.</div>';
  }
  function chartColors(){const dark=document.documentElement.dataset.theme==='dark';return {paper_bgcolor:'rgba(0,0,0,0)',plot_bgcolor:'rgba(0,0,0,0)',
    font:{color:dark?'#e9f2ff':'#243356',family:'Inter, Segoe UI, sans-serif',size:11},
    xaxis:{gridcolor:dark?'#35445a':'#d5dfeb',zeroline:false,linecolor:dark?'#667690':'#8593a9'},
    yaxis:{gridcolor:dark?'#35445a':'#d5dfeb',zeroline:false,linecolor:dark?'#667690':'#8593a9'}};
  }
  const colors=['#06cbd2','#a54ce8','#efa14a','#4fbb8a','#ef6494','#347fe6'];
  function drawLine(target, data, title) {
    if(target==='plot' && $('plot').querySelector('img')) $('plot').replaceChildren();
    const palette=chartColors();
    const traces=data.traces.map((t,i)=>({x:t.x||data.x,y:t.y,type:'scatter',mode:'lines',name:t.name,
      line:{width:3,color:colors[i%colors.length]},hovertemplate:'%{x:.3f} · %{y:.5g}<extra>%{fullData.name}</extra>'}));
    const layout={...palette,xaxis:{...palette.xaxis,title:{text:data.xlabel||''},automargin:true},
      yaxis:{...palette.yaxis,title:{text:data.ylabel||''},automargin:true},
      margin:{l:72,r:12,t:12,b:56},legend:{orientation:'h',y:1.13,x:0.03,bgcolor:'rgba(0,0,0,0)',font:{size:11}},
      showlegend:true,hovermode:'x unified'};
    Plotly.react(target,traces,layout,{responsive:true,displaylogo:false,scrollZoom:true});
    if(title)text('chart-title',title);
  }
  function plotUrl(uid, metric, download=false) {
    const name = document.documentElement.dataset.theme==='dark' ? 'dark' : 'light';
    return `/api/runs/${encodeURIComponent(uid)}/plot/${encodeURIComponent(metric)}.png?theme=${name}${download?'&download=1':''}`;
  }
  function showSavedImage(container, uid, metric, caption) {
    const holder=$(container);
    if (typeof Plotly!=='undefined' && holder.classList.contains('js-plotly-plot')) Plotly.purge(holder);
    const img=document.createElement('img');
    img.className='generated-plot-image';
    img.alt=caption;
    img.loading='eager';
    img.decoding='async';
    img.addEventListener('error',()=>{
      if(img.isConnected) holder.textContent='Unable to load the generated plot. Check the Runner diagnostics log or try Interactive mode.';
    });
    img.src=plotUrl(uid,metric);
    holder.replaceChildren(img);
  }
  function updatePlotModeButton() {
    const button=$('chart-interactive');
    button.setAttribute('aria-pressed',state.plotMode==='interactive'?'true':'false');
    button.textContent=state.plotMode==='interactive'?'Saved PNG':'Interactive';
    button.title=state.plotMode==='interactive'?'Show saved Matplotlib image':'Switch to interactive Plotly chart';
  }
  async function renderMetric() {
    if(state.view!=='plots'||state.section!=='simulation') return;
    const ticket=++state.plotRequest;
    const metric=state.metric;
    $$('#metric-tabs button').forEach(el=>el.classList.toggle('selected',el.dataset.metric===metric));
    updatePlotModeButton();
    const uid=state.active;
    if(!uid) {
      if(metric==='kernels'||metric==='atmosphere'){
        if(typeof Plotly!=='undefined')Plotly.purge('plot');
        $('plot').textContent='Run a simulation to visualize these results.';
        return;
      }
      try {
        const demo=await fetch('/static/sample_nadir.json').then(r=>r.json());
        if(ticket!==state.plotRequest) return;
        if(metric==='transmittance'){demo.traces=demo.transmittance;demo.ylabel='Transmittance (–)';}
        drawLine('plot',demo,metric==='transmittance'?'Sample Nadir · Transmittance':'Sample Nadir · Brightness Temperature');
      }catch(e){if(ticket===state.plotRequest)$('plot').textContent='Sample chart unavailable.';}
      text('chart-note','Illustrative preview · Run a simulation for real plots.');
      return;
    }
    let run;
    try {
      // Do not use state.runs here: it can contain a queued/running snapshot even after completion.
      run=await api(`/api/runs/${uid}`);
    } catch(e) {if(ticket===state.plotRequest)$('plot').textContent=e.message;return;}
    if(ticket!==state.plotRequest||uid!==state.active||metric!==state.metric) return;
    if(run.status!=='completed') {
      if(typeof Plotly!=='undefined')Plotly.purge('plot');
      $('plot').textContent=run.status==='failed'?'Simulation failed: '+(run.error||'See diagnostics.'):
        run.status==='interrupted'?'Simulation interrupted.':
        'Simulation '+run.status+' ('+(run.stage||'starting')+'). Results will appear automatically.';
      return;
    }
    const title=metric==='kernels'?'Temperature Sensitivity · Kernel':
      metric==='atmosphere'?'Atmospheric Profiles':
      metric==='transmittance'?'Transmittance · '+run.geometry:
      (run.parameters.output==='bbt'?'Brightness Temperature':'Radiance')+' · '+run.geometry;
    text('chart-title',title);
    if(metric==='kernels'&&!run.parameters.kernel){
      if(typeof Plotly!=='undefined')Plotly.purge('plot');
      $('plot').textContent='This run did not include kernel calculations. Enable “Calculate Kernels” for a future run.';
      text('chart-note','No kernel.tab for this run.');return;
    }
    if(state.plotMode==='image') {
      showSavedImage('plot',uid,metric,title);
      text('chart-note',`Saved Matplotlib plot · run ${uid} · ${run.warnings} warning(s) · Use Interactive to zoom and inspect values.`);
      return;
    }
    try {
      if(metric==='kernels') {
        const info=await api(`/api/runs/${uid}/kernel`);
        if(ticket!==state.plotRequest||uid!==state.active||metric!==state.metric) return;
        const pal=chartColors();
        if($('plot').querySelector('img'))$('plot').replaceChildren();
        Plotly.react('plot',[{type:'heatmap',x:info.x,y:info.y,z:info.z,colorscale:'Teal',colorbar:{title:'Sensitivity'}}],
          {...pal,xaxis:{...pal.xaxis,title:info.xlabel},yaxis:{...pal.yaxis,title:info.ylabel},margin:{l:65,b:48,t:10,r:10}},
          {responsive:true,displaylogo:false});
      } else if(metric==='atmosphere') {
        const data=await api(`/api/runs/${uid}/chart/atmosphere`);
        if(ticket!==state.plotRequest||uid!==state.active||metric!==state.metric) return;
        drawAtmosphericChart('plot',data);
      } else {
        const data=await api(`/api/runs/${uid}/chart/${metric}`);
        if(ticket!==state.plotRequest||uid!==state.active||metric!==state.metric) return;
        drawLine('plot',data,title);
      }
      text('chart-note',`Interactive JURASSIC results · run ${uid} · ${run.warnings} warning(s)`);
    } catch(e) {
      // A Plotly failure should never prevent the user from seeing the saved scientific plot.
      if(ticket===state.plotRequest){
        state.plotMode='image';updatePlotModeButton();
        showSavedImage('plot',uid,metric,title);
        toast('Interactive chart unavailable; displaying saved PNG instead.',true);
      }
    }
  }
  function drawAtmosphericChart(id, data) {
    if(id==='plot' && $('plot').querySelector('img')) $('plot').replaceChildren();
    const palette=chartColors();
    const lines=data.traces.map((t,i)=>({type:'scatter',mode:'lines',x:t.x,y:t.y,name:t.name,line:{width:2.5,color:colors[i%colors.length]}}));
    Plotly.react(id,lines,{...palette,xaxis:{...palette.xaxis,title:'Atmospheric value (mixed units)'},yaxis:{...palette.yaxis,title:'Altitude (km)'},
      margin:{l:57,r:12,t:24,b:53},legend:{orientation:'h',y:1.15,font:{size:10}}},{responsive:true,displaylogo:false});
  }
  async function renderAtmosphere(){
    if(state.view!=='atmosphere') return;
    if(!state.active){text('atmosphere-plot','Complete a simulation to view its atmospheric profile.');return;}
    try {
      const run=await api(`/api/runs/${state.active}`);
      if(state.view!=='atmosphere'||state.active!==run.id)return;
      if(run.status!=='completed'){
        text('atmosphere-plot','Atmospheric plots will appear after the simulation completes.');return;
      }
      showSavedImage('atmosphere-plot',run.id,'atmosphere','Temperature, pressure and gas profiles');
    }catch(e){text('atmosphere-plot',e.message);}
  }
  async function renderTable(){
    if(state.view!=='table')return;
    const r=state.runs.find(x=>x.id===state.active);
    if(!r||r.status!=='completed') return;
    try{
      const data=await api(`/api/runs/${state.active}/table`);
      $('results-table').innerHTML='<table><thead><tr>'+data.headers.map((v,i)=>`<th title="${escapeHtml(v)}">${escapeHtml('C'+(i+1))}</th>`).join('')+
        '</tr></thead><tbody>'+data.rows.map(row=>'<tr>'+row.map(v=>`<td>${escapeHtml(Number(v).toPrecision(6))}</td>`).join('')+'</tr>').join('')+'</tbody></table>'+
        `<p class="subtle">Showing ${data.preview_rows} rows · Hover a column heading for the full variable name. Download CSV for all rows.</p>`;
    }catch(e){text('results-table',e.message);}
  }
  async function renderConfig() {
    const r=state.runs.find(x=>x.id===state.active);
    if(!r){if(!$('config-text').textContent.trim())text('config-text','Create or preview a configuration.');return;}
    try{
      const response=await fetch(`/api/runs/${r.id}/download/config.ctl`);
      if(response.ok)text('config-text',await response.text());
    }catch(e){/* temporary worker startup */}
  }
  async function renderLogs() {
    const r=state.runs.find(x=>x.id===state.active);
    if(!r){text('log-text','No run selected.');return;}
    try{const data=await api(`/api/runs/${r.id}/log?step=${encodeURIComponent($('log-step').value)}`);text('log-text',data.text);}
    catch(e){text('log-text',e.message);}
  }
  function renderCatalog(){
    const rows=[];
    Object.entries(state.catalog).forEach(([id,ds])=>{
      Object.entries(ds.channels).forEach(([channel,info])=>{
        Object.entries(info.gases).forEach(([gas,entry])=>{
          rows.push(`<tr><td>${escapeHtml(ds.label)}</td><td>${escapeHtml(channel)}</td><td>${escapeHtml(gas)}</td>
            <td>${entry.data?'<span class="badge">Numeric table</span>':'<span class="badge">Header only</span>'}</td>
            <td>${info.filter?'Yes':'Missing'}</td><td>${(entry.bytes/1048576).toFixed(2)} MB</td></tr>`);
        });
      });
    });
    $('catalog-table').innerHTML=rows.length?'<table><thead><tr><th>Dataset</th><th>Channel (cm⁻¹)</th><th>Gas</th><th>Table</th><th>Filter</th><th>Size</th></tr></thead><tbody>'+rows.join('')+'</tbody></table>':'<div class="empty-state">No compatible lookup tables discovered. Check the upstream installation.</div>';
  }
  async function uploadDataset() {
    const btn=$('upload-dataset-button');btn.disabled=true;
    try{
      const label=$('dataset-label').value.trim(),file=$('dataset-upload').files[0];
      if(!file)throw new Error('Select a ZIP dataset bundle');
      const form=new FormData();form.append('label',label);form.append('dataset',file);
      const out=await api('/api/catalog/upload','POST',form);
      text('dataset-upload-message',`Registered ${out.label} (${out.channels} channels).`);
      await loadCatalog();renderCatalog();toast('Dataset registered. Check gas/channel compatibility before running.');
    }catch(e){toast(e.message,true);}finally{btn.disabled=false;}
  }
  async function doCompare() {
    const ids=$$('.compare-choice:checked').map(e=>e.value);
    if(ids.length<2||ids.length>4){toast('Select 2–4 completed runs.',true);return;}
    try{
      const info=await api('/api/compare?ids='+ids.join(',')+'&metric='+$('compare-metric').value);
      drawLine('compare-plot',info);
      $('compare-plot').dataset.ready='true';
      text('compare-note',`Overlay of ${ids.length} completed JURASSIC simulations. Compare physical units and geometry before interpreting differences.`);
    }catch(e){toast(e.message,true);}
  }
  async function refreshSystem() {
    try{
      const r=await api('/api/system');
      $('sidebar-status').textContent=r.ready?'JURASSIC ready':'Build required';
      document.querySelector('.status-dot').classList.toggle('ready',r.ready);
      $('system-overview').innerHTML=`<p>Git revision: <code>${escapeHtml(r.revision)}</code> · Registered datasets: ${r.datasets} · Queued/running: ${r.queued}</p>`+
        Object.entries(r.executables).map(([name,good])=>`<span class="system-binary ${good?'ok':''}">${good?'✓':'✗'} ${escapeHtml(name)}</span>`).join('');
    }catch(e){text('system-overview',e.message);text('sidebar-status','Offline');}
  }
  async function loadCatalog(){
    const result=await api('/api/catalog');state.catalog=result;
    if(!state.catalog[state.dataset]){state.dataset=Object.keys(state.catalog)[0]||'';resetSelections();}
    paintSetup();
  }
  async function exportCsv(){
    if(!state.active){toast('Run a simulation first.',true);return;}
    const r=state.runs.find(x=>x.id===state.active);
    if(r?.status==='completed')window.location.href=`/api/runs/${state.active}/download/rad.csv`;
  }
  function bind() {
    $('light-btn').onclick=()=>theme('light');$('dark-btn').onclick=()=>theme('dark');
    $$('.nav-item, .nav-submenu button').forEach(el=>el.onclick=()=>section(el.dataset.section));
    $$('.flow-step').forEach(el=>el.onclick=()=>flow(el.dataset.flow));
    $$('#geometry-options button').forEach(el=>el.onclick=()=>setGeometry(el.dataset.geometry));
    $('dataset-select').onchange=()=>{state.dataset=$('dataset-select').value;resetSelections();paintSetup();};
    $('basic-mode').onclick=()=>{$('advanced-fields').hidden=true;$('basic-mode').classList.add('selected');$('advanced-mode').classList.remove('selected');};
    $('advanced-mode').onclick=()=>{$('advanced-fields').hidden=false;$('advanced-mode').classList.add('selected');$('basic-mode').classList.remove('selected');};
    $('add-channel-shortcut').onclick=()=>section('datasets');
    $('run-button').onclick=runSimulation;
    $('preview-btn').onclick=async()=>{try{await previewConfig();setView('config');toast('Configuration preview generated.');}catch(e){toast(e.message,true);}};
    $('save-config-btn').onclick=async()=>{try{const v=await previewConfig();saveFile(v.text,`${state.geometry}.ctl`);}catch(e){toast(e.message,true);}};
    $('config-studio-preview').onclick=async()=>{try{await previewConfig();}catch(e){toast(e.message,true);}};
    $('config-studio-save').onclick=async()=>{try{const r=await previewConfig();saveFile(r.text,`${state.geometry}.ctl`);}catch(e){toast(e.message,true);}};
    $('config-file').onchange=async()=>{
      const file=$('config-file').files[0];if(!file)return;
      const form=new FormData();form.append('configuration',file);
      try{
        const imported=await api('/api/config/import','POST',form);
        const v=imported.values;
        if(v.gases?.length)state.gases=v.gases.filter(g=>allGases(currentDataset()).includes(g));
        if(v.channels?.length)state.channels=v.channels.map(Number).map(n=>n.toFixed(4)).filter(c=>currentDataset().channels[c]);
        if(v.output){state.output=v.output;$$('input[name="output"]').forEach(el=>el.checked=el.value===v.output);}
        if(v.rayds)$('rayds').value=v.rayds;if(v.raydz)$('raydz').value=v.raydz;
        if(v.refrac)$('refrac').checked=v.refrac!=='0';
        paintSetup();toast('Supported values imported. Review the selected dataset and gas/channel compatibility.');
      }catch(e){toast(e.message,true);}
    };
    $$('#content-tabs button').forEach(el=>el.onclick=()=>setView(el.dataset.view));
    $$('#metric-tabs button').forEach(el=>el.onclick=()=>{state.metric=el.dataset.metric;renderMetric();});
    $('log-step').onchange=renderLogs;
    $('copy-ctl').onclick=()=>navigator.clipboard?.writeText($('config-text').textContent).then(()=>toast('Configuration copied.')).catch(()=>toast('Clipboard access unavailable.',true));
    $('chart-interactive').onclick=()=>{state.plotMode=state.plotMode==='interactive'?'image':'interactive';renderMetric();};
    $('chart-expand').onclick=()=>{const el=$('plot');if(el.requestFullscreen)el.requestFullscreen().catch(()=>toast('Fullscreen unavailable.',true));};
    $('chart-download').onclick=()=>{
      if(state.active){
        window.location.href=plotUrl(state.active,state.metric,true);
      }else if(typeof Plotly!=='undefined'){
        Plotly.downloadImage('plot',{format:'png',filename:'jurassic-preview',width:1200,height:700});
      }else toast('Complete a simulation to download a plot.',true);
    };
    $('csv-download').onclick=exportCsv;
    $('upload-dataset-button').onclick=uploadDataset;
    $('compare-button').onclick=doCompare;
    $('history-refresh').onclick=()=>refreshHistory().catch(e=>toast(e.message,true));
    $$('input[name="output"]').forEach(el=>el.addEventListener('change',()=>{state.output=el.value;}));
  }
  async function init(){
    bind();theme(localStorage.getItem('jurassic-theme')==='dark'?'dark':'light');
    try {await loadCatalog();await refreshSystem();await refreshHistory();}
    catch(e){toast(e.message,true);}
    if(Object.keys(state.catalog).length){setGeometry('nadir');}else{text('compatibility-note','No tables found. Check the upstream tests/data directory.');}
    await renderMetric();
  }
  document.addEventListener('DOMContentLoaded',init);
})();
