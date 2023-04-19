function setupEmail()
%% Email setup
    from_address = 'breuerflume@gmail.com';
    pass = 'lfvxzaqmulqfaumi'; % email password is flumepass
    setpref('Internet', 'E_mail', from_address);
    setpref('Internet', 'SMTP_Username', from_address);
    setpref('Internet', 'SMTP_Password', pass);
    setpref('Internet', 'SMTP_Server',   'smtp.gmail.com');
    props = java.lang.System.getProperties;
    props.setProperty('mail.smtp.auth',                'true');  % Note: 'true' as a string, not a logical value!
    props.setProperty('mail.smtp.starttls.enable',     'true');  % Note: 'true' as a string, not a logical value!
    props.setProperty('mail.smtp.socketFactory.port',  '465');   % Note: '465'  as a string, not a numeric value!
    props.setProperty('mail.smtp.socketFactory.class', 'javax.net.ssl.SSLSocketFactory');
end