var ze=Object.create;var U=Object.defineProperty;var je=Object.getOwnPropertyDescriptor;var Be=Object.getOwnPropertyNames;var We=Object.getPrototypeOf,Ve=Object.prototype.hasOwnProperty;var k=(e=>typeof require<"u"?require:typeof Proxy<"u"?new Proxy(e,{get:(t,r)=>(typeof require<"u"?require:t)[r]}):e)(function(e){if(typeof require<"u")return require.apply(this,arguments);throw Error('Dynamic require of "'+e+'" is not supported')});var qe=(e,t)=>()=>(t||e((t={exports:{}}).exports,t),t.exports),Ye=(e,t)=>{for(var r in t)U(e,r,{get:t[r],enumerable:!0})},Ge=(e,t,r,o)=>{if(t&&typeof t=="object"||typeof t=="function")for(let a of Be(t))!Ve.call(e,a)&&a!==r&&U(e,a,{get:()=>t[a],enumerable:!(o=je(t,a))||o.enumerable});return e};var Je=(e,t,r)=>(r=e!=null?ze(We(e)):{},Ge(t||!e||!e.__esModule?U(r,"default",{value:e,enumerable:!0}):r,e));var ae=qe(()=>{});var M={};Ye(M,{createEndpoint:()=>N,expose:()=>E,finalizer:()=>T,proxy:()=>I,proxyMarker:()=>j,releaseProxy:()=>ee,transfer:()=>oe,transferHandlers:()=>v,windowEndpoint:()=>nt,wrap:()=>A});var j=Symbol("Comlink.proxy"),N=Symbol("Comlink.endpoint"),ee=Symbol("Comlink.releaseProxy"),T=Symbol("Comlink.finalizer"),C=Symbol("Comlink.thrown"),te=e=>typeof e=="object"&&e!==null||typeof e=="function",Xe={canHandle:e=>te(e)&&e[j],serialize(e){let{port1:t,port2:r}=new MessageChannel;return E(e,t),[r,[r]]},deserialize(e){return e.start(),A(e)}},Ke={canHandle:e=>te(e)&&C in e,serialize({value:e}){let t;return e instanceof Error?t={isError:!0,value:{message:e.message,name:e.name,stack:e.stack}}:t={isError:!1,value:e},[t,[]]},deserialize(e){throw e.isError?Object.assign(new Error(e.value.message),e.value):e.value}},v=new Map([["proxy",Xe],["throw",Ke]]);function Qe(e,t){for(let r of e)if(t===r||r==="*"||r instanceof RegExp&&r.test(t))return!0;return!1}function E(e,t=globalThis,r=["*"]){t.addEventListener("message",function o(a){if(!a||!a.data)return;if(!Qe(r,a.origin)){console.warn(`Invalid origin '${a.origin}' for comlink proxy`);return}let{id:i,type:n,path:c}=Object.assign({path:[]},a.data),s=(a.data.argumentList||[]).map(x),u;try{let l=c.slice(0,-1).reduce((m,y)=>m[y],e),f=c.reduce((m,y)=>m[y],e);switch(n){case"GET":u=f;break;case"SET":l[c.slice(-1)[0]]=x(a.data.value),u=!0;break;case"APPLY":u=f.apply(l,s);break;case"CONSTRUCT":{let m=new f(...s);u=I(m)}break;case"ENDPOINT":{let{port1:m,port2:y}=new MessageChannel;E(e,y),u=oe(m,[m])}break;case"RELEASE":u=void 0;break;default:return}}catch(l){u={value:l,[C]:0}}Promise.resolve(u).catch(l=>({value:l,[C]:0})).then(l=>{let[f,m]=R(l);t.postMessage(Object.assign(Object.assign({},f),{id:i}),m),n==="RELEASE"&&(t.removeEventListener("message",o),re(t),T in e&&typeof e[T]=="function"&&e[T]())}).catch(l=>{let[f,m]=R({value:new TypeError("Unserializable return value"),[C]:0});t.postMessage(Object.assign(Object.assign({},f),{id:i}),m)})}),t.start&&t.start()}function Ze(e){return e.constructor.name==="MessagePort"}function re(e){Ze(e)&&e.close()}function A(e,t){return z(e,[],t)}function F(e){if(e)throw new Error("Proxy has been released and is not useable")}function ne(e){return P(e,{type:"RELEASE"}).then(()=>{re(e)})}var L=new WeakMap,_="FinalizationRegistry"in globalThis&&new FinalizationRegistry(e=>{let t=(L.get(e)||0)-1;L.set(e,t),t===0&&ne(e)});function et(e,t){let r=(L.get(t)||0)+1;L.set(t,r),_&&_.register(e,t,e)}function tt(e){_&&_.unregister(e)}function z(e,t=[],r=function(){}){let o=!1,a=new Proxy(r,{get(i,n){if(F(o),n===ee)return()=>{tt(a),ne(e),o=!0};if(n==="then"){if(t.length===0)return{then:()=>a};let c=P(e,{type:"GET",path:t.map(s=>s.toString())}).then(x);return c.then.bind(c)}return z(e,[...t,n])},set(i,n,c){F(o);let[s,u]=R(c);return P(e,{type:"SET",path:[...t,n].map(l=>l.toString()),value:s},u).then(x)},apply(i,n,c){F(o);let s=t[t.length-1];if(s===N)return P(e,{type:"ENDPOINT"}).then(x);if(s==="bind")return z(e,t.slice(0,-1));let[u,l]=Z(c);return P(e,{type:"APPLY",path:t.map(f=>f.toString()),argumentList:u},l).then(x)},construct(i,n){F(o);let[c,s]=Z(n);return P(e,{type:"CONSTRUCT",path:t.map(u=>u.toString()),argumentList:c},s).then(x)}});return et(a,e),a}function rt(e){return Array.prototype.concat.apply([],e)}function Z(e){let t=e.map(R);return[t.map(r=>r[0]),rt(t.map(r=>r[1]))]}var ie=new WeakMap;function oe(e,t){return ie.set(e,t),e}function I(e){return Object.assign(e,{[j]:!0})}function nt(e,t=globalThis,r="*"){return{postMessage:(o,a)=>e.postMessage(o,r,a),addEventListener:t.addEventListener.bind(t),removeEventListener:t.removeEventListener.bind(t)}}function R(e){for(let[t,r]of v)if(r.canHandle(e)){let[o,a]=r.serialize(e);return[{type:"HANDLER",name:t,value:o},a]}return[{type:"RAW",value:e},ie.get(e)||[]]}function x(e){switch(e.type){case"HANDLER":return v.get(e.name).deserialize(e.value);case"RAW":return e.value}}function P(e,t,r){return new Promise(o=>{let a=it();e.addEventListener("message",function i(n){!n.data||!n.data.id||n.data.id!==a||(e.removeEventListener("message",i),o(n.data))}),e.start&&e.start(),e.postMessage(Object.assign({id:a},t),r)})}function it(){return new Array(4).fill(0).map(()=>Math.floor(Math.random()*Number.MAX_SAFE_INTEGER).toString(16)).join("-")}var ot=Object.create,V=Object.defineProperty,at=Object.getOwnPropertyDescriptor,st=Object.getOwnPropertyNames,lt=Object.getPrototypeOf,ct=Object.prototype.hasOwnProperty,d=(e,t)=>V(e,"name",{value:t,configurable:!0}),ce=(e=>typeof k<"u"?k:typeof Proxy<"u"?new Proxy(e,{get:(t,r)=>(typeof k<"u"?k:t)[r]}):e)(function(e){if(typeof k<"u")return k.apply(this,arguments);throw new Error('Dynamic require of "'+e+'" is not supported')}),ue=(e,t)=>()=>(t||e((t={exports:{}}).exports,t),t.exports),ut=(e,t,r,o)=>{if(t&&typeof t=="object"||typeof t=="function")for(let a of st(t))!ct.call(e,a)&&a!==r&&V(e,a,{get:()=>t[a],enumerable:!(o=at(t,a))||o.enumerable});return e},ft=(e,t,r)=>(r=e!=null?ot(lt(e)):{},ut(t||!e||!e.__esModule?V(r,"default",{value:e,enumerable:!0}):r,e)),mt=ue((e,t)=>{(function(r,o){"use strict";typeof define=="function"&&define.amd?define("stackframe",[],o):typeof e=="object"?t.exports=o():r.StackFrame=o()})(e,function(){"use strict";function r(p){return!isNaN(parseFloat(p))&&isFinite(p)}d(r,"_isNumber");function o(p){return p.charAt(0).toUpperCase()+p.substring(1)}d(o,"_capitalize");function a(p){return function(){return this[p]}}d(a,"_getter");var i=["isConstructor","isEval","isNative","isToplevel"],n=["columnNumber","lineNumber"],c=["fileName","functionName","source"],s=["args"],u=["evalOrigin"],l=i.concat(n,c,s,u);function f(p){if(p)for(var g=0;g<l.length;g++)p[l[g]]!==void 0&&this["set"+o(l[g])](p[l[g]])}d(f,"StackFrame"),f.prototype={getArgs:function(){return this.args},setArgs:function(p){if(Object.prototype.toString.call(p)!=="[object Array]")throw new TypeError("Args must be an Array");this.args=p},getEvalOrigin:function(){return this.evalOrigin},setEvalOrigin:function(p){if(p instanceof f)this.evalOrigin=p;else if(p instanceof Object)this.evalOrigin=new f(p);else throw new TypeError("Eval Origin must be an Object or StackFrame")},toString:function(){var p=this.getFileName()||"",g=this.getLineNumber()||"",b=this.getColumnNumber()||"",O=this.getFunctionName()||"";return this.getIsEval()?p?"[eval] ("+p+":"+g+":"+b+")":"[eval]:"+g+":"+b:O?O+" ("+p+":"+g+":"+b+")":p+":"+g+":"+b}},f.fromString=d(function(p){var g=p.indexOf("("),b=p.lastIndexOf(")"),O=p.substring(0,g),De=p.substring(g+1,b).split(","),Q=p.substring(b+1);if(Q.indexOf("@")===0)var H=/@(.+?)(?::(\d+))?(?::(\d+))?$/.exec(Q,""),$e=H[1],He=H[2],Ue=H[3];return new f({functionName:O,args:De||void 0,fileName:$e,lineNumber:He||void 0,columnNumber:Ue||void 0})},"StackFrame$$fromString");for(var m=0;m<i.length;m++)f.prototype["get"+o(i[m])]=a(i[m]),f.prototype["set"+o(i[m])]=function(p){return function(g){this[p]=!!g}}(i[m]);for(var y=0;y<n.length;y++)f.prototype["get"+o(n[y])]=a(n[y]),f.prototype["set"+o(n[y])]=function(p){return function(g){if(!r(g))throw new TypeError(p+" must be a Number");this[p]=Number(g)}}(n[y]);for(var w=0;w<c.length;w++)f.prototype["get"+o(c[w])]=a(c[w]),f.prototype["set"+o(c[w])]=function(p){return function(g){this[p]=String(g)}}(c[w]);return f})}),pt=ue((e,t)=>{(function(r,o){"use strict";typeof define=="function"&&define.amd?define("error-stack-parser",["stackframe"],o):typeof e=="object"?t.exports=o(mt()):r.ErrorStackParser=o(r.StackFrame)})(e,d(function(r){"use strict";var o=/(^|@)\S+:\d+/,a=/^\s*at .*(\S+:\d+|\(native\))/m,i=/^(eval@)?(\[native code])?$/;return{parse:d(function(n){if(typeof n.stacktrace<"u"||typeof n["opera#sourceloc"]<"u")return this.parseOpera(n);if(n.stack&&n.stack.match(a))return this.parseV8OrIE(n);if(n.stack)return this.parseFFOrSafari(n);throw new Error("Cannot parse given Error object")},"ErrorStackParser$$parse"),extractLocation:d(function(n){if(n.indexOf(":")===-1)return[n];var c=/(.+?)(?::(\d+))?(?::(\d+))?$/,s=c.exec(n.replace(/[()]/g,""));return[s[1],s[2]||void 0,s[3]||void 0]},"ErrorStackParser$$extractLocation"),parseV8OrIE:d(function(n){var c=n.stack.split(`
`).filter(function(s){return!!s.match(a)},this);return c.map(function(s){s.indexOf("(eval ")>-1&&(s=s.replace(/eval code/g,"eval").replace(/(\(eval at [^()]*)|(,.*$)/g,""));var u=s.replace(/^\s+/,"").replace(/\(eval code/g,"(").replace(/^.*?\s+/,""),l=u.match(/ (\(.+\)$)/);u=l?u.replace(l[0],""):u;var f=this.extractLocation(l?l[1]:u),m=l&&u||void 0,y=["eval","<anonymous>"].indexOf(f[0])>-1?void 0:f[0];return new r({functionName:m,fileName:y,lineNumber:f[1],columnNumber:f[2],source:s})},this)},"ErrorStackParser$$parseV8OrIE"),parseFFOrSafari:d(function(n){var c=n.stack.split(`
`).filter(function(s){return!s.match(i)},this);return c.map(function(s){if(s.indexOf(" > eval")>-1&&(s=s.replace(/ line (\d+)(?: > eval line \d+)* > eval:\d+:\d+/g,":$1")),s.indexOf("@")===-1&&s.indexOf(":")===-1)return new r({functionName:s});var u=/((.*".+"[^@]*)?[^@]*)(?:@)/,l=s.match(u),f=l&&l[1]?l[1]:void 0,m=this.extractLocation(s.replace(u,""));return new r({functionName:f,fileName:m[0],lineNumber:m[1],columnNumber:m[2],source:s})},this)},"ErrorStackParser$$parseFFOrSafari"),parseOpera:d(function(n){return!n.stacktrace||n.message.indexOf(`
`)>-1&&n.message.split(`
`).length>n.stacktrace.split(`
`).length?this.parseOpera9(n):n.stack?this.parseOpera11(n):this.parseOpera10(n)},"ErrorStackParser$$parseOpera"),parseOpera9:d(function(n){for(var c=/Line (\d+).*script (?:in )?(\S+)/i,s=n.message.split(`
`),u=[],l=2,f=s.length;l<f;l+=2){var m=c.exec(s[l]);m&&u.push(new r({fileName:m[2],lineNumber:m[1],source:s[l]}))}return u},"ErrorStackParser$$parseOpera9"),parseOpera10:d(function(n){for(var c=/Line (\d+).*script (?:in )?(\S+)(?:: In function (\S+))?$/i,s=n.stacktrace.split(`
`),u=[],l=0,f=s.length;l<f;l+=2){var m=c.exec(s[l]);m&&u.push(new r({functionName:m[3]||void 0,fileName:m[2],lineNumber:m[1],source:s[l]}))}return u},"ErrorStackParser$$parseOpera10"),parseOpera11:d(function(n){var c=n.stack.split(`
`).filter(function(s){return!!s.match(o)&&!s.match(/^Error created at/)},this);return c.map(function(s){var u=s.split("@"),l=this.extractLocation(u.pop()),f=u.shift()||"",m=f.replace(/<anonymous function(: (\w+))?>/,"$2").replace(/\([^)]*\)/g,"")||void 0,y;f.match(/\(([^)]*)\)/)&&(y=f.replace(/^[^(]+\(([^)]*)\)$/,"$1"));var w=y===void 0||y==="[arguments not available]"?void 0:y.split(",");return new r({functionName:m,args:w,fileName:l[0],lineNumber:l[1],columnNumber:l[2],source:s})},this)},"ErrorStackParser$$parseOpera11")}},"ErrorStackParser"))}),dt=ft(pt()),h=typeof process=="object"&&typeof process.versions=="object"&&typeof process.versions.node=="string"&&typeof process.browser>"u",fe=h&&typeof module<"u"&&typeof module.exports<"u"&&typeof ce<"u"&&typeof __dirname<"u",yt=h&&!fe,gt=typeof Deno<"u",me=!h&&!gt,wt=me&&typeof window=="object"&&typeof document=="object"&&typeof document.createElement=="function"&&typeof sessionStorage=="object"&&typeof importScripts!="function",ht=me&&typeof importScripts=="function"&&typeof self=="object",Tt=typeof navigator=="object"&&typeof navigator.userAgent=="string"&&navigator.userAgent.indexOf("Chrome")==-1&&navigator.userAgent.indexOf("Safari")>-1,pe,B,de,se,q;async function Y(){if(!h||(pe=(await import("node:url")).default,se=await import("node:fs"),q=await import("node:fs/promises"),de=(await import("node:vm")).default,B=await import("node:path"),G=B.sep,typeof ce<"u"))return;let e=se,t=await import("node:crypto"),r=await Promise.resolve().then(()=>Je(ae(),1)),o=await import("node:child_process"),a={fs:e,crypto:t,ws:r,child_process:o};globalThis.require=function(i){return a[i]}}d(Y,"initNodeModules");function ye(e,t){return B.resolve(t||".",e)}d(ye,"node_resolvePath");function ge(e,t){return t===void 0&&(t=location),new URL(e,t).toString()}d(ge,"browser_resolvePath");var W;h?W=ye:W=ge;var G;h||(G="/");function we(e,t){return e.startsWith("file://")&&(e=e.slice(7)),e.includes("://")?{response:fetch(e)}:{binary:q.readFile(e).then(r=>new Uint8Array(r.buffer,r.byteOffset,r.byteLength))}}d(we,"node_getBinaryResponse");function he(e,t){let r=new URL(e,location);return{response:fetch(r,t?{integrity:t}:{})}}d(he,"browser_getBinaryResponse");var $;h?$=we:$=he;async function ve(e,t){let{response:r,binary:o}=$(e,t);if(o)return o;let a=await r;if(!a.ok)throw new Error(`Failed to load '${e}': request failed.`);return new Uint8Array(await a.arrayBuffer())}d(ve,"loadBinaryFile");var D;if(wt)D=d(async e=>await import(e),"loadScript");else if(ht)D=d(async e=>{try{globalThis.importScripts(e)}catch(t){if(t instanceof TypeError)await import(e);else throw t}},"loadScript");else if(h)D=be;else throw new Error("Cannot determine runtime environment");async function be(e){e.startsWith("file://")&&(e=e.slice(7)),e.includes("://")?de.runInThisContext(await(await fetch(e)).text()):await import(pe.pathToFileURL(e).href)}d(be,"nodeLoadScript");async function ke(e){if(h){await Y();let t=await q.readFile(e,{encoding:"utf8"});return JSON.parse(t)}else return await(await fetch(e)).json()}d(ke,"loadLockFile");async function xe(){if(fe)return __dirname;let e;try{throw new Error}catch(o){e=o}let t=dt.default.parse(e)[0].fileName;if(yt){let o=await import("node:path");return(await import("node:url")).fileURLToPath(o.dirname(t))}let r=t.lastIndexOf(G);if(r===-1)throw new Error("Could not extract indexURL path from pyodide module location");return t.slice(0,r)}d(xe,"calculateDirname");function Ee(e){let t=e.FS,r=e.FS.filesystems.MEMFS,o=e.PATH,a={DIR_MODE:16895,FILE_MODE:33279,mount:function(i){if(!i.opts.fileSystemHandle)throw new Error("opts.fileSystemHandle is required");return r.mount.apply(null,arguments)},syncfs:async(i,n,c)=>{try{let s=a.getLocalSet(i),u=await a.getRemoteSet(i),l=n?u:s,f=n?s:u;await a.reconcile(i,l,f),c(null)}catch(s){c(s)}},getLocalSet:i=>{let n=Object.create(null);function c(l){return l!=="."&&l!==".."}d(c,"isRealDir");function s(l){return f=>o.join2(l,f)}d(s,"toAbsolute");let u=t.readdir(i.mountpoint).filter(c).map(s(i.mountpoint));for(;u.length;){let l=u.pop(),f=t.stat(l);t.isDir(f.mode)&&u.push.apply(u,t.readdir(l).filter(c).map(s(l))),n[l]={timestamp:f.mtime,mode:f.mode}}return{type:"local",entries:n}},getRemoteSet:async i=>{let n=Object.create(null),c=await vt(i.opts.fileSystemHandle);for(let[s,u]of c)s!=="."&&(n[o.join2(i.mountpoint,s)]={timestamp:u.kind==="file"?(await u.getFile()).lastModifiedDate:new Date,mode:u.kind==="file"?a.FILE_MODE:a.DIR_MODE});return{type:"remote",entries:n,handles:c}},loadLocalEntry:i=>{let n=t.lookupPath(i).node,c=t.stat(i);if(t.isDir(c.mode))return{timestamp:c.mtime,mode:c.mode};if(t.isFile(c.mode))return n.contents=r.getFileDataAsTypedArray(n),{timestamp:c.mtime,mode:c.mode,contents:n.contents};throw new Error("node type not supported")},storeLocalEntry:(i,n)=>{if(t.isDir(n.mode))t.mkdirTree(i,n.mode);else if(t.isFile(n.mode))t.writeFile(i,n.contents,{canOwn:!0});else throw new Error("node type not supported");t.chmod(i,n.mode),t.utime(i,n.timestamp,n.timestamp)},removeLocalEntry:i=>{var n=t.stat(i);t.isDir(n.mode)?t.rmdir(i):t.isFile(n.mode)&&t.unlink(i)},loadRemoteEntry:async i=>{if(i.kind==="file"){let n=await i.getFile();return{contents:new Uint8Array(await n.arrayBuffer()),mode:a.FILE_MODE,timestamp:n.lastModifiedDate}}else{if(i.kind==="directory")return{mode:a.DIR_MODE,timestamp:new Date};throw new Error("unknown kind: "+i.kind)}},storeRemoteEntry:async(i,n,c)=>{let s=i.get(o.dirname(n)),u=t.isFile(c.mode)?await s.getFileHandle(o.basename(n),{create:!0}):await s.getDirectoryHandle(o.basename(n),{create:!0});if(u.kind==="file"){let l=await u.createWritable();await l.write(c.contents),await l.close()}i.set(n,u)},removeRemoteEntry:async(i,n)=>{await i.get(o.dirname(n)).removeEntry(o.basename(n)),i.delete(n)},reconcile:async(i,n,c)=>{let s=0,u=[];Object.keys(n.entries).forEach(function(m){let y=n.entries[m],w=c.entries[m];(!w||t.isFile(y.mode)&&y.timestamp.getTime()>w.timestamp.getTime())&&(u.push(m),s++)}),u.sort();let l=[];if(Object.keys(c.entries).forEach(function(m){n.entries[m]||(l.push(m),s++)}),l.sort().reverse(),!s)return;let f=n.type==="remote"?n.handles:c.handles;for(let m of u){let y=o.normalize(m.replace(i.mountpoint,"/")).substring(1);if(c.type==="local"){let w=f.get(y),p=await a.loadRemoteEntry(w);a.storeLocalEntry(m,p)}else{let w=a.loadLocalEntry(m);await a.storeRemoteEntry(f,y,w)}}for(let m of l)if(c.type==="local")a.removeLocalEntry(m);else{let y=o.normalize(m.replace(i.mountpoint,"/")).substring(1);await a.removeRemoteEntry(f,y)}}};e.FS.filesystems.NATIVEFS_ASYNC=a}d(Ee,"initializeNativeFS");var vt=d(async e=>{let t=[];async function r(a){for await(let i of a.values())t.push(i),i.kind==="directory"&&await r(i)}d(r,"collect"),await r(e);let o=new Map;o.set(".",e);for(let a of t){let i=(await e.resolve(a)).join("/");o.set(i,a)}return o},"getFsHandles");function Pe(e){let t={noImageDecoding:!0,noAudioDecoding:!0,noWasmDecoding:!1,preRun:Ce(e),quit(r,o){throw t.exited={status:r,toThrow:o},o},print:e.stdout,printErr:e.stderr,arguments:e.args,API:{config:e},locateFile:r=>e.indexURL+r,instantiateWasm:Le(e.indexURL)};return t}d(Pe,"createSettings");function Se(e){return function(t){let r="/";try{t.FS.mkdirTree(e)}catch(o){console.error(`Error occurred while making a home directory '${e}':`),console.error(o),console.error(`Using '${r}' for a home directory instead`),e=r}t.FS.chdir(e)}}d(Se,"createHomeDirectory");function Oe(e){return function(t){Object.assign(t.ENV,e)}}d(Oe,"setEnvironment");function Fe(e){return t=>{for(let r of e)t.FS.mkdirTree(r),t.FS.mount(t.FS.filesystems.NODEFS,{root:r},r)}}d(Fe,"mountLocalDirectories");function Te(e){let t=ve(e);return r=>{let o=r._py_version_major(),a=r._py_version_minor();r.FS.mkdirTree("/lib"),r.FS.mkdirTree(`/lib/python${o}.${a}/site-packages`),r.addRunDependency("install-stdlib"),t.then(i=>{r.FS.writeFile(`/lib/python${o}${a}.zip`,i)}).catch(i=>{console.error("Error occurred while installing the standard library:"),console.error(i)}).finally(()=>{r.removeRunDependency("install-stdlib")})}}d(Te,"installStdlib");function Ce(e){let t;return e.stdLibURL!=null?t=e.stdLibURL:t=e.indexURL+"python_stdlib.zip",[Te(t),Se(e.env.HOME),Oe(e.env),Fe(e._node_mounts),Ee]}d(Ce,"getFileSystemInitializationFuncs");function Le(e){let{binary:t,response:r}=$(e+"pyodide.asm.wasm");return function(o,a){return async function(){try{let i;r?i=await WebAssembly.instantiateStreaming(r,o):i=await WebAssembly.instantiate(await t,o);let{instance:n,module:c}=i;typeof WasmOffsetConverter<"u"&&(wasmOffsetConverter=new WasmOffsetConverter(wasmBinary,c)),a(n,c)}catch(i){console.warn("wasm instantiation failed!"),console.warn(i)}}(),{}}}d(Le,"getInstantiateWasmFunc");var le="0.26.1";async function J(e={}){await Y();let t=e.indexURL||await xe();t=W(t),t.endsWith("/")||(t+="/"),e.indexURL=t;let r={fullStdLib:!1,jsglobals:globalThis,stdin:globalThis.prompt?globalThis.prompt:void 0,lockFileURL:t+"pyodide-lock.json",args:[],_node_mounts:[],env:{},packageCacheDir:t,packages:[],enableRunUntilComplete:!1},o=Object.assign(r,e);o.env.HOME||(o.env.HOME="/home/pyodide");let a=Pe(o),i=a.API;if(i.lockFilePromise=ke(o.lockFileURL),typeof _createPyodideModule!="function"){let l=`${o.indexURL}pyodide.asm.js`;await D(l)}let n;if(e._loadSnapshot){let l=await e._loadSnapshot;ArrayBuffer.isView(l)?n=l:n=new Uint8Array(l),a.noInitialRun=!0,a.INITIAL_MEMORY=n.length}let c=await _createPyodideModule(a);if(a.exited)throw a.exited.toThrow;if(e.pyproxyToStringRepr&&i.setPyProxyToStringMethod(!0),i.version!==le)throw new Error(`Pyodide version does not match: '${le}' <==> '${i.version}'. If you updated the Pyodide version, make sure you also updated the 'indexURL' parameter passed to loadPyodide.`);c.locateFile=l=>{throw new Error("Didn't expect to load any more file_packager files!")};let s;n&&(s=i.restoreSnapshot(n));let u=i.finalizeBootstrap(s);return i.sys.path.insert(0,i.config.env.HOME),u.version.includes("dev")||i.setCdnUrl(`https://cdn.jsdelivr.net/pyodide/v${u.version}/full/`),i._pyodide.set_excepthook(),await i.packageIndexReady,i.initializeStreams(o.stdin,o.stdout,o.stderr),u}d(J,"loadPyodide");function X(e){return typeof ImageBitmap<"u"&&e instanceof ImageBitmap}function S(e,t,r,...o){return e==null||X(e)||e instanceof ArrayBuffer||ArrayBuffer.isView(e)?e:t(e)?r(e,...o):Array.isArray(e)?e.map(a=>S(a,t,r,...o)):typeof e=="object"?Object.fromEntries(Object.entries(e).map(([a,i])=>[a,S(i,t,r,...o)])):e}function bt(e){return e&&e[Symbol.toStringTag]=="PyProxy"}function _e(e){return e&&!!e[N]}function kt(e){return e&&typeof e=="object"&&"_comlinkProxy"in e&&"ptr"in e}function xt(e){return e&&e[Symbol.toStringTag]=="Map"}function K(e){if(_e(e))return!0;if(e==null||e instanceof ArrayBuffer||ArrayBuffer.isView(e))return!1;if(e instanceof Array)return e.some(t=>K(t));if(typeof e=="object")return Object.entries(e).some(([t,r])=>K(r))}var Re={},Ne={canHandle:bt,serialize(e){let t=self.pyodide._module.PyProxy_getPtr(e);Re[t]=e;let{port1:r,port2:o}=new MessageChannel;return E(e,r),[[o,t],[o]]},deserialize([e,t]){e.start();let r=A(e);return new Proxy(r,{get:(a,i)=>i==="_ptr"?t:a[i]})}},Ae={canHandle:K,serialize(e){return[S(e,_e,t=>({_comlinkProxy:!0,ptr:t._ptr})),[]]},deserialize(e){return S(e,kt,t=>Re[t.ptr])}},Ie={canHandle:X,serialize(e){if(e.width==0&&e.height==0){let t=new OffscreenCanvas(1,1);t.getContext("2d"),e=t.transferToImageBitmap()}return[e,[e]]},deserialize(e){return e}},Me={canHandle:xt,serialize(e){return[Object.fromEntries(e.entries()),[]]},deserialize(e){return e}};var Et={mkdir(e){self.pyodide._FS.mkdir(e)},writeFile(e,t){self.pyodide._FS.writeFile(e,t)}};async function Pt(e){return self.pyodide=await J(e),self.pyodide.registerComlink(M),self.pyodide._FS=self.pyodide.FS,self.pyodide.FS=Et,v.set("PyProxy",Ne),v.set("Comlink",Ae),v.set("ImageBitmap",Ie),v.set("Map",Me),I(self.pyodide)}E({init:Pt});
/*! Bundled license information:

comlink/dist/esm/comlink.mjs:
  (**
   * @license
   * Copyright 2019 Google LLC
   * SPDX-License-Identifier: Apache-2.0
   *)
*/