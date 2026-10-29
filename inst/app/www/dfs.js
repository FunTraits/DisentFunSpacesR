// DisentFunSpaces application: navigation from the home page and step marks
$(document).on('click keydown', '[data-go]', function (e) {
  if (e.type === 'keydown' && e.key !== 'Enter' && e.key !== ' ') return;
  e.preventDefault();
  Shiny.setInputValue('go_tab', String($(this).data('go')), {priority: 'event'});
});
$(document).on('shiny:connected', function () {
  Shiny.addCustomMessageHandler('dfs_steps', function (msg) {
    var done = [].concat(msg.done || []);
    $('.step-n[data-step]').each(function () {
      var s = Number($(this).data('step'));
      $(this).toggleClass('done', done.indexOf(s) >= 0);
    });
  });
});
